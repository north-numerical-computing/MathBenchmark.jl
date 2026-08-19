using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using MathBenchmark
MathBenchmark.run_mathbenchmark()
