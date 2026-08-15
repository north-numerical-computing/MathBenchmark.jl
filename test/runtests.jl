using Test
using MathBenchmark
const Err = MathBenchmark.Err
const Config = MathBenchmark.Config

function with_config(f, contents::AbstractString)
    mktempdir() do dir
        file = joinpath(dir, "config.json")
        write(file, contents)
        f(file)
    end
end

@testset "MathBenchmark" begin
    @testset "config" begin
        valid = """
        {
            "second" : { "format" : "binary64", "rounding" : "RN", "fastmath" : 0, "search" : "hours" },
            "first"  : { "format" : "binary16", "rounding" : "RZ", "fastmath" : 1, "search" : "exhaustive" },
            "third"  : { "format" : "binary32", "rounding" : "RU", "fastmath" : false, "search" : 12345 }
        }
        """
        tasks = with_config(Config.read_config, valid)
        # Order of the file is preserved.
        @test [t.name for t in tasks] == ["second", "first", "third"]
        @test tasks[1].format == "binary64" && tasks[1].search == "hours" &&
              tasks[1].rounding == "RN" && tasks[1].fastmath === false
        @test tasks[2].format == "binary16" && tasks[2].search == "exhaustive" &&
              tasks[2].rounding == "RZ" && tasks[2].fastmath === true
        @test tasks[3].search === 12345 && tasks[3].fastmath === false

        task(fields) = """{ "t" : { $fields } }"""
        good = """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0, "search" : "seconds" """
        @test with_config(Config.read_config, task(good))[1].name == "t"
        # Every invalid or missing field is an error, not a silently different run.
        for bad in (
                """ "format" : "binary128", "rounding" : "RN", "fastmath" : 0, "search" : "seconds" """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0, "search" : "hour" """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0, "search" : 0 """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0, "search" : -5 """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0, "search" : true """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 0 """,
                """ "format" : "binary32", "rounding" : "RX", "fastmath" : 0, "search" : "seconds" """,
                """ "format" : "binary32", "fastmath" : 0, "search" : "seconds" """,
                """ "format" : "binary32", "rounding" : "RN", "fastmath" : 2, "search" : "seconds" """,
                """ "format" : "binary32", "rounding" : "RN", "search" : "seconds" """,
                """ "rounding" : "RN", "fastmath" : 0, "search" : "seconds" """,
            )
            @test_throws ArgumentError with_config(Config.read_config, task(bad))
        end
        @test_throws ArgumentError with_config(Config.read_config, """{ "t" : 3 }""")
        @test_throws ArgumentError with_config(Config.read_config, """[1, 2]""")
    end

    @testset "resolve_function" begin
        @test Err.resolve_function("exp", false) === exp
        @test Err.resolve_function("exp", true) === Base.FastMath.exp_fast
        @test Err.resolve_function("sin", true) === Base.FastMath.sin_fast
        # No fast variant: the plain function is used.
        @test Err.resolve_function("sinpi", true) === sinpi
        @test Err.resolve_function("tanpi", true) === tanpi
        @test_throws ArgumentError Err.resolve_function("nosuchfunction", false)
        # The reference is always the plain BigFloat implementation.
        @test Err.reference_function("exp") === exp
        # The fast variants actually differ from the plain functions somewhere.
        @test any(Base.FastMath.exp_fast(x) != exp(x) for x in range(-700.0, 700.0, length=1001))
    end

    @testset "ordinals ($T)" for T in (Float16, Float32, Float64)
        specials = T[-Inf, -floatmax(T), -1.0, -floatmin(T), -nextfloat(zero(T)), -0.0,
                     0.0, nextfloat(zero(T)), floatmin(T), 1.0, floatmax(T), Inf]
        for x in specials
            @test Err.float_from_ordinal(T, Err.ordinal(x)) === x
        end
        @test Err.ordinal(T(-0.0)) == -1
        @test Err.ordinal(T(0.0)) == 0
        @test Err.ordinal(-nextfloat(zero(T))) == -2
        # Consecutive floats have consecutive ordinals.
        for x in T[-1.0, -floatmin(T), -nextfloat(zero(T)), floatmin(T), 1.0, floatmax(T)]
            @test Err.ordinal(nextfloat(x)) == Err.ordinal(x) + 1
            @test Err.ordinal(prevfloat(x)) == Err.ordinal(x) - 1
        end
        # Unlike nextfloat/prevfloat, ordinals do not skip one of the zeros.
        @test Err.ordinal(nextfloat(T(-0.0))) == 1
        @test Err.ordinal(prevfloat(T(0.0))) == -2
        @test Err.number_of_floats(T(1), T(1)) == 1
        @test Err.number_of_floats(T(-0.0), T(0.0)) == 2
        @test Err.number_of_floats(-floatmax(T), floatmax(T)) ==
              2 * (Int128(reinterpret(Base.uinttype(T), floatmax(T))) + 1)
        # Number of floats in an interval by brute force (Float16 only, it is small).
        if T === Float16
            lo, hi = T(-2.5), T(3.25)
            n = 1  # nextfloat skips +0.0 when coming from -0.0
            x = lo
            while x <= hi
                n += 1
                x = nextfloat(x)
            end
            @test Err.number_of_floats(lo, hi) == n
        end
    end

    @testset "split_range" begin
        for (lo, hi, nchunks) in ((0, 99, 7), (-50, 49, 8), (3, 3, 4), (-5, 5, 100),
                                  (Err.ordinal(-floatmax(Float64)), Err.ordinal(floatmax(Float64)), 128))
            chunks = MathBenchmark.split_range(lo, hi, nchunks)
            @test length(chunks) == min(nchunks, Int128(hi) - Int128(lo) + 1)
            @test first(chunks)[1] == lo && last(chunks)[2] == hi
            @test all(a <= b for (a, b) in chunks)
            @test all(chunks[i][2] + 1 == chunks[i + 1][1] for i in 1:length(chunks) - 1)
            lens = [b - a + 1 for (a, b) in chunks]
            @test maximum(lens) - minimum(lens) <= 1
        end
        @test_throws ArgumentError MathBenchmark.split_range(5, 4, 3)
    end

    @testset "ulp error" begin
        setprecision(BigFloat, 127) do
            # y is exactly the reference: no error.
            @test Err.get_ulp_error(1.5, big"1.5") == 0
            # Half an ulp above 1 (ulp = 2^-52 in [1, 2)).
            @test Err.get_ulp_error(1.0, big"1.0" + big"2.0"^-53) == 0.5
            # Half an ulp below 1: the ulp is that of [0.5, 1) even though the
            # reference rounds up to the power of two.
            @test Err.get_ulp_error(1.0, big"1.0" - big"2.0"^-54) == 0.5
            @test Err.get_ulp_error(prevfloat(1.0), big"1.0" - big"2.0"^-54) == 0.5
            # Exact power of two: ulp of the binade above.
            @test Err.get_ulp_error(nextfloat(1.0), big"1.0") == 1
            # Subnormals: constant spacing, also just below floatmin.
            @test Err.get_ulp_error(nextfloat(0.0), big"1.5" * BigFloat(nextfloat(0.0))) == 0.5
            @test Err.get_ulp_error(floatmin(Float64), BigFloat(floatmin(Float64)) - BigFloat(nextfloat(0.0)) / 2) == 0.5
            # Zero reference.
            @test Err.get_ulp_error(nextfloat(0.0), big"0.0") == 1
            # Sign and other formats.
            @test Err.get_ulp_error(-1.0f0, -big"1.0" - big"2.0"^-25) == 0.25
            @test Err.get_ulp_error(Float16(2.0), big"2.0" + big"2.0"^-9 * 3) == 3
        end
    end

    @testset "kernels" begin
        T = Float16
        f = sin
        ref = Err.BigFloatReference(sin)
        lo, hi = T(-3.0), T(3.0)
        setprecision(BigFloat, Err.reference_precision(T)) do
            # Exhaustive over a small interval, in one task and in many, against a
            # direct computation with the oracle.
            r1 = MathBenchmark.test_domain(ref, f, T, lo, hi, 1; nchunks=1)
            r2 = MathBenchmark.test_domain(ref, f, T, lo, hi, 1; nchunks=13)
            errs = [Float64(Err.get_ulp_error(f(x), sin(BigFloat(x))), RoundToZero)
                    for x in (Err.float_from_ordinal(T, k) for k in Err.ordinal(lo):Err.ordinal(hi))]
            @test r1.ntests == r2.ntests == length(errs) == Err.number_of_floats(lo, hi)
            @test r1.max_error == r2.max_error == maximum(errs)
            @test r1.input === r2.input && r1.output === r2.output == f(r1.input)
            @test r1.reference == r2.reference == sin(BigFloat(r1.input))
            @test r1.ninfs == r2.ninfs == 0
            # Fixed step: the sampled set is ordinal(lo) + j*step regardless of chunking.
            step = 7
            r3 = MathBenchmark.test_domain(ref, f, T, lo, hi, step; nchunks=1)
            r4 = MathBenchmark.test_domain(ref, f, T, lo, hi, step; nchunks=5)
            sampled = Err.ordinal(lo):step:Err.ordinal(hi)
            @test r3.ntests == r4.ntests == length(sampled)
            @test r3.max_error == r4.max_error == maximum(errs[1:step:end])
            @test r3.input === r4.input
            # Special inputs and infinities.
            r5 = Err.test_values(ref, f, T[0.5, 1.5, 2.5])
            @test r5.ntests == 3 && r5.max_error == maximum(errs[Err.ordinal(T(x)) - Err.ordinal(lo) + 1] for x in (0.5, 1.5, 2.5))
            r6 = Err.test_values(Err.BigFloatReference(exp), exp, T[1.0, 12.0, 20.0])
            @test r6.ntests == 3 && r6.ninfs == 2 && r6.input === T(1.0)
            # Combining results keeps the worst case and sums the counters.
            c = Err.combine(r5, r6)
            @test c.ntests == 6 && c.ninfs == 2 && c.max_error == max(r5.max_error, r6.max_error)
            @test Err.combine(r6, Err.Result{T}()).max_error == r6.max_error
        end
    end

    @testset "run_mathbenchmark" begin
        config = """
        {
            "b16" : { "format" : "binary16", "rounding" : "RN", "fastmath" : 0, "search" : 2 },
            "b64" : { "format" : "binary64", "rounding" : "RN", "fastmath" : 1, "search" : 2 }
        }
        """
        with_config(config) do file
            outdir = joinpath(dirname(file), "out")
            MathBenchmark.run_mathbenchmark(file; output_dir=outdir)
            for name in ("b16", "b64"), prefix in ("", "HEX_")
                lines = readlines(joinpath(outdir, "$(prefix)$name.txt"))
                @test length(lines) == 1 + length(MathBenchmark.Functions.functions_dict["binary64"])
                @test startswith(lines[1], "Function")
                rows = split.(lines[2:end])
                @test all(0 <= parse(Float64, r[2]) < 10 for r in rows)
                @test allunique(r[1] for r in rows)
            end
        end
    end
end
