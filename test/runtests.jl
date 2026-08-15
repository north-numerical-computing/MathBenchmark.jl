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
end
