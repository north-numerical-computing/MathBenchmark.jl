module MathBenchmark

include("file_io.jl")
using .Config

include("functions.jl")
using .Functions

include("error_calculator.jl")
using .Err

using Printf

# Number of tasks spawned per thread when testing a function. More chunks than
# threads let the scheduler balance the load: the cost of the reference
# computation varies with the magnitude of the input.
const CHUNKS_PER_THREAD = 8

# Number of tests per thread used to estimate the speed of a function, and to
# warm up (compile) the kernels before timing them.
const CALIBRATION_TESTS = 100_000
const WARMUP_TESTS = 100

# Time budget, in nanoseconds, of the time-based search strategies.
const TIME_BUDGET_NS = Dict("seconds" => 10^9,
                            "minutes" => 60 * 10^9,
                            "hours"   => 3600 * 10^9,
                            "days"    => 24 * 3600 * 10^9)

"""
    split_range(k_lo, k_hi, nchunks) -> Vector{Tuple{Int128, Int128}}

Split the integer range `k_lo:k_hi` into at most `nchunks` contiguous, non-empty
chunks of nearly equal length.
"""
function split_range(k_lo::Integer, k_hi::Integer, nchunks::Integer)
    lo = Int128(k_lo)
    n = Int128(k_hi) - lo + 1
    n >= 1 || throw(ArgumentError("empty range $k_lo:$k_hi"))
    nchunks = clamp(Int128(nchunks), 1, n)
    return [(lo + div((i - 1) * n, nchunks), lo + div(i * n, nchunks) - 1) for i in 1:nchunks]
end

"""
    test_domain(reference, f, ::Type{T}, lo, hi, step; nchunks) -> Result{T}

Test `f` on the values `x` of `[lo, hi]` whose ordinal is `ordinal(lo) + j * step`,
i.e. every `step`-th value starting from `lo` (`step = 1` is exhaustive). The
domain is split into `nchunks` chunks tested by concurrent tasks whose results
are combined.
"""
function test_domain(reference, f, ::Type{T}, lo::T, hi::T, step::Integer;
                     nchunks::Integer = CHUNKS_PER_THREAD * Threads.nthreads()) where {T}
    k_lo = Err.ordinal(lo)
    step = Int128(step)
    # With many threads, Julia (seen on 1.11, 1.12 and nightly with 96 threads)
    # can deadlock when a garbage collection is requested while a thread holds
    # MPFR's internal cache lock (taken e.g. to compute pi for acos): the lock
    # holder waits for the collection at the safepoint inside
    # jl_gc_counted_malloc, which MPFR uses for its temporaries, while the
    # threads waiting for the lock sit in a ccall and never reach a safepoint,
    # so the collection cannot start. The kernels do not allocate Julia
    # objects, so the collector is simply disabled while the tasks run.
    gc_enabled = GC.enable(false)
    try
        tasks = map(split_range(k_lo, Err.ordinal(hi), nchunks)) do (a, b)
            # First sampled ordinal in this chunk.
            first = k_lo + cld(a - k_lo, step) * step
            Threads.@spawn Err.test_range(reference, f, T, first, b, step)
        end
        return reduce(Err.combine, fetch.(tasks))
    finally
        GC.enable(gc_enabled)
    end
end

"""
    tests_per_thread_in_budget(reference, f, ::Type{T}, lo, hi, budget_ns;
                               ntests = CALIBRATION_TESTS) -> Int

Estimate how many tests per thread of `f` fit in `budget_ns` nanoseconds by timing
`ntests` tests per thread (after a short warm-up run, so that the compilation of
the kernels for `f` is not counted).
"""
function tests_per_thread_in_budget(reference, f, ::Type{T}, lo::T, hi::T, budget_ns;
                                    ntests::Integer = CALIBRATION_TESTS) where {T}
    nthreads = Threads.nthreads()
    num_floats = Err.number_of_floats(lo, hi)
    stride(n) = max(cld(num_floats, n), 1)
    test_domain(reference, f, T, lo, hi, stride(WARMUP_TESTS * nthreads))
    elapsed = @elapsed test_domain(reference, f, T, lo, hi, stride(ntests * nthreads))
    return max(floor(Int, budget_ns / (elapsed * 1e9 / ntests)), 1)
end

"""
    run_function(func_name, ::Type{T}, lo, hi, search, fastmath_on) -> Result{T}

Test the function `func_name` on its input domain `[lo, hi]` with the search
strategy `search` (see `Config.BenchmarkTask`), then on its special inputs.
"""
function run_function(func_name::AbstractString, ::Type{T}, lo::T, hi::T, search, fastmath_on::Bool) where {T}
    f = Err.resolve_function(func_name, fastmath_on)
    reference = Err.MPFRReference(func_name, Err.reference_precision(T))
    nthreads = Threads.nthreads()
    num_floats = Err.number_of_floats(lo, hi)

    # Number of tests to run in the domain, and the corresponding stride. A
    # budget allowing more tests than there are values falls back to the
    # exhaustive search.
    ntests = if search == "exhaustive"
        num_floats
    elseif search isa Integer
        # An integer is the number of tests per thread.
        Int128(search) * nthreads
    else
        Int128(tests_per_thread_in_budget(reference, f, T, lo, hi, TIME_BUDGET_NS[search])) * nthreads
    end
    exhaustive = ntests >= num_floats
    step = exhaustive ? Int128(1) : cld(num_floats, ntests)
    ntests = cld(num_floats, step)

    if exhaustive
        @printf("Running %d tests (search strategy exhaustive) for the function \
                 %s with %d threads \n", ntests, func_name, nthreads)
    else
        @printf("Running %d tests (search strategy \"%s\") for the function \
                 %s with %d threads \n", ntests, search, func_name, nthreads)
    end
    flush(stdout)

    result = test_domain(reference, f, T, lo, hi, step)

    # Inputs known to be hard for other libraries.
    special_inputs = Functions.special_inputs(T, func_name)
    if !isempty(special_inputs)
        result = Err.combine(result, Err.test_values(reference, f, special_inputs))
    end
    return result
end

"""
    run_mathbenchmark(config_file = "config.json"; output_dir = "output")

Run the tasks described in the JSON file `config_file` and write one pair of result
files (`<task>.txt` and `HEX_<task>.txt`) per task in `output_dir`.
"""
function run_mathbenchmark(config_file::AbstractString = "config.json";
                           output_dir::AbstractString = "output")

    # Read and validate the testing tasks specified in the json file. Invalid
    # tasks throw an error before anything is run.
    tasks = Config.read_config(config_file)
    mkpath(output_dir)

    for task in tasks
        run_task(task, output_dir)
    end
end

function run_task(task::Config.BenchmarkTask, output_dir::AbstractString)
    task_name = task.name
    T = Err.float_type(task.format)
    printstyled("Running task: $task_name\n", color=:blue)
    printstyled("Format: $(task.format) Search: $(task.search) Rounding: \
                $(task.rounding) Fastmath: $(task.fastmath)\n", color=:green)

    # Results table formatting. Each task specified in the JSON file has a
    # result .txt file named accordingly, with the input and output printed
    # with 17 significant digits (enough to round-trip a binary64 value) and
    # the reference with 21, and a HEX_ file with the bits of input and output.
    file = joinpath(output_dir, "$task_name.txt")
    file_hex = joinpath(output_dir, "HEX_$task_name.txt")
    write(file, Printf.format(Printf.Format("%-10s %15s %30s %30s %30s %20s %20s %20s\n"),
                              "Function", "ULPs", "Input", "Output", "MPFR", "Tests", "Infs", "Failures"))
    write(file_hex, Printf.format(Printf.Format("%-10s %15s %30s %30s %20s %20s %20s\n"),
                                  "Function", "ULPs", "Input", "Output", "Tests", "Infs", "Failures"))
    fe = Printf.Format("%-10s %15.10f %30.16e %30.16e %30.20e %20d %20d %20d\n")
    fe_hex = Printf.Format("%-10s %15.10f %#30x %#30x %20d %20d %20d\n")

    # Loop through the functions list of a particular format, in alphabetical
    # order so that result files are comparable across runs and Julia versions.
    domains = Functions.functions_dict[task.format]
    for func_name in sort!(collect(keys(domains)))
        lo, hi = domains[func_name]
        r = run_function(func_name, T, lo, hi, task.search, task.fastmath)

        # Report the maximum error and the corresponding values to the output files.
        max_error = trunc(r.max_error, digits=10)
        open(file, "a") do io
            Printf.format(io, fe, func_name, max_error, Float64(r.input), Float64(r.output),
                          r.reference, r.ntests, r.ninfs, r.nfailures)
        end
        open(file_hex, "a") do io
            Printf.format(io, fe_hex, func_name, max_error,
                          reinterpret(Base.uinttype(T), r.input),
                          reinterpret(Base.uinttype(T), r.output),
                          r.ntests, r.ninfs, r.nfailures)
        end
    end
end

end
