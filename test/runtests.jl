using Test
using MathBenchmark
const Err = MathBenchmark.Err

@testset "MathBenchmark" begin
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
        @test any(Base.FastMath.exp_fast(x) != exp(x) for x in range(-700.0, 700.0, length=100_001))
    end
end
