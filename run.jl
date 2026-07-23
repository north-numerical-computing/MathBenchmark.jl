using Logging
Logging.disable_logging(Logging.Warn)

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using MathBenchmark
MathBenchmark.run_benchmark()

