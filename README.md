# A Test Bench for Measuring the Errors of Mathematical Functions of Julia

This repository contains the data produced in [Sec. 4.3, 1] and the testing code used to generate the results.

In order to run Julia mathematical function accuracy tests start by modifying [config.json](./config.json). In this file a configuration of the formats, rounding modes, and input space search strategies is provided to the test bench.

Then, run the script [MathBenchmark.jl](./MathBenchmark.jl) to start the testing specified in [config.json](./config.json). The code is threaded and can be run by

```
julia --threads=auto MathBenchmark.jl
```

Once the testing completes, the results can be found in the [output](./output) directory which will contain the text files for each of the test cases specified in [config.jl](./config.jl). The obtained results can be compared with the accuracy of thirteen mathematical function libraries that were analysed by Gladman, Innocente, Mather, and Zimmermann [2].

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

The `search` field can be set to `seconds`, `minutes`, `hours`, `days`, or `exhaustive`; the first four options determine the approximate duration of testing one mathematical function. For example, the `seconds` option will mean that the input domain of each function is traversed in a fixed-size step which allows the testing of each function to take approximately one second.
Options `exhaustive` tells the benchmark to test all possible inputs in the function's input domain, provided exactly in [functions.jl](./functions.jl) for a particular floating-point format.

### References

 [1] M. Mikaitis and T. Rizyal, [*Accuracy of Mathematical Functions in Julia*](https://arxiv.org/pdf/2509.05666). 	arXiv:2509.05666 [cs.MS]. Sep., 2025.

 [2] B. Gladman, V. Innocente, J. Mather, and P. Zimmermann. [*Accuracy of mathematical functions in single, double, double extended, and quadruple precision*](https://members.loria.fr/PZimmermann/papers/accuracy.pdf). Preprint. Aug., 2025.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENCE](./LICENCE)).
