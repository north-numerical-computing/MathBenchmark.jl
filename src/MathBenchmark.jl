module MathBenchmark

include("file_io.jl")
using .Config

include("functions.jl")
using .Functions

include("error_calculator.jl")
using .Err

using Printf
using Format
using Logging

# Number of tasks spawned per thread when testing a function. More chunks than
# threads let the scheduler balance the load: the cost of the reference
# computation varies with the magnitude of the input.
const CHUNKS_PER_THREAD = 8

# Number of tests per thread used to estimate the speed of a function.
const CALIBRATION_TESTS = 100_000

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
    tasks = map(split_range(k_lo, Err.ordinal(hi), nchunks)) do (a, b)
        # First sampled ordinal in this chunk.
        first = k_lo + cld(a - k_lo, step) * step
        Threads.@spawn Err.test_range(reference, f, T, first, b, step)
    end
    return reduce(Err.combine, fetch.(tasks))
end

"""
    tests_per_thread_in_budget(reference, f, ::Type{T}, lo, hi, budget_ns) -> Int

Estimate how many tests per thread of `f` fit in `budget_ns` nanoseconds by timing
`CALIBRATION_TESTS` tests per thread.
"""
function tests_per_thread_in_budget(reference, f, ::Type{T}, lo::T, hi::T, budget_ns) where {T}
    nthreads = Threads.nthreads()
    ntests = CALIBRATION_TESTS * nthreads
    step = max(cld(Err.number_of_floats(lo, hi), ntests), 1)
    elapsed = @elapsed test_domain(reference, f, T, lo, hi, step)
    return floor(Int, budget_ns / (elapsed * 1e9 / CALIBRATION_TESTS))
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

    # Results table formatting. Each task specied in the JSON file has
    # a result .txt file named accordingly.
    fe = FormatExpr("{1:<10s} {2:>15s} {3:>30s} {4:>30s} {5:>30s} {6:>20s} {7:>20s}\n")
    result_table_head = format(fe, "Function", "ULPs", "Input", "Output",
                               "MPFR", "Tests", "Infs")
    fe_hex = FormatExpr("{1:<10s} {2:>15s} {3:>30s} {4:>30s} {5:>20s}\n")
    result_table_head_hex = format(fe_hex, "Function", "ULPs", "Input", "Output",
                                   "Tests")
    file = joinpath(output_dir, "$task_name.txt")
    file_hex = joinpath(output_dir, "HEX_$task_name.txt")
    write(file, result_table_head)
    write(file_hex, result_table_head_hex)
    fe = FormatExpr("{1:<10s} {2:>15.10f} {3:>30.15e} \
                    {4:>30.15e} {5:>30.15e} {6:>20d} {7:>20d}\n")
    fe_hex = FormatExpr("{1:<10s} {2:>15.10f} {3:>#30x} {4:>#30x} \
                    {5:>20d}\n")

    # Loop through the functions list of a particular format.
    for (func_name, (lo, hi)) in Functions.functions_dict[task.format]
        r = run_function(func_name, T, lo, hi, task.search, task.fastmath)

        # Report the maximum error and the corresponding values to the output files.
        max_error = trunc(r.max_error, digits=10)
        line = format(fe, func_name, max_error, Float64(r.input), Float64(r.output),
                      Float64(r.reference), r.ntests, r.ninfs)
        open(file, "a") do io
            write(io, line)
        end
        line = format(fe_hex, func_name, max_error,
                      reinterpret(Base.uinttype(T), r.input),
                      reinterpret(Base.uinttype(T), r.output),
                      r.ntests)
        open(file_hex, "a") do io
            write(io, line)
        end
    end
end

end
