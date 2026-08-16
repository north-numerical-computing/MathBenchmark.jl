# A Test Bench for Measuring the Errors of Mathematical Functions of Julia

This repository contains the data produced in [Sec. 4.3, 1] and the testing code used to generate the results.

In order to run Julia mathematical function accuracy tests start by modifying [config.json](./config.json). In this file a configuration of the formats, rounding modes, and input space search strategies is provided to the test bench.

Then, run the script [run.jl](./run.jl) to start the testing specified in [config.json](./config.json). The code is threaded and can be run by

```
julia --threads=auto run.jl
```

The work is dominated by MPFR, so using one thread per physical core (rather than per hardware thread) is usually best; on machines with many cores, enabling the concurrent GC sweeper (`--gcthreads=N,1`) also helps. From Julia, the same can be done with `MathBenchmark.run_mathbenchmark(config_file; output_dir)`.

Once the testing completes, the results can be found in the [output](./output) directory which will contain the text files for each of the test cases specified in [config.json](./config.json): `<task>.txt` with, for each function, the maximum error in ulps found, the input for which it occurs, the output of Julia and the correctly rounded MPFR reference (`HEX_<task>.txt` gives the bits of input and output), the number of tests, the number of infinities returned by the function (`Infs`) and the number of results that could not be compared to the reference (`Failures`: `Inf`/`NaN` where the correctly rounded result is finite, a finite result where it is infinite, or an infinity of the wrong sign). The obtained results can be compared with the accuracy of thirteen mathematical function libraries that were analysed by Gladman, Innocente, Mather, and Zimmermann [2].

Here is an example configuration JSON file:

```
{
    "test-binary16RN-exhaustive-nofastmath" : {
        "format" : "binary16",
        "rounding" : "RN",
        "fastmath" : 0,
        "search" : "exhaustive"
    },
    "test-binary32RN-seconds-nofastmath" : {
        "format" : "binary32",
        "rounding" : "RN",
        "fastmath" : 0,
        "search" : "hours"
    },
    "test-binary64RN-seconds-nofastmath" : {
        "format" : "binary64",
        "rounding" : "RN",
        "fastmath" : 0,
        "search" : "hours"
    }
}
```

The `format` field can be set to `binary16`, `binary32`, or `binary64`. The `rounding` field is a feature that may be available in the future; this can be set to `RN`, `RZ`, `RU`, or `RD`, but at present it will not take effect because Julia does not provide mathematical functions with separate rounding modes.

The field `fastmath` can be set either to `0` or `1` to turn the Julia's fastmath feature off or on, respectively.

The `search` field can be set to `seconds`, `minutes`, `hours`, `days`, `exhaustive`, or to an integer value; the first four options determine the approximate duration of testing one mathematical function. For example, the `seconds` option will mean that the input domain of each function is traversed in a fixed-size step which allows the testing of each function to take approximately one second. This is a very rough approximation - it is achieved by measuring the time to test the particular function on 100000 inputs per thread and extrapolating. The `search` field also accepts an integer value for the number of tests to undertake per thread: the domain is then sampled with a fixed stride so that the total number of tests is that value times the number of threads. In both cases the domain is traversed exhaustively if it contains fewer values than the number of tests.
The search setting `exhaustive` tells the benchmark to test all possible inputs in the function's input domain, provided exactly in [functions.jl](./src/functions.jl) for a particular floating-point format.

The tests can be run with `julia --project -e 'using Pkg; Pkg.test()'`.

### References

 [1] M. Mikaitis and T. Rizyal, [*Accuracy of Mathematical Functions in Julia*](https://arxiv.org/pdf/2509.05666). 	arXiv:2509.05666 [cs.MS]. Sep., 2025.

 [2] B. Gladman, V. Innocente, J. Mather, and P. Zimmermann. [*Accuracy of mathematical functions in single, double, double extended, and quadruple precision*](https://members.loria.fr/PZimmermann/papers/accuracy.pdf). Preprint. Aug., 2025.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENCE](./LICENCE)).
