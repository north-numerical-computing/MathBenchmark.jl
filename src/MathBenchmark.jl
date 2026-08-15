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


mutable struct PaddedFloat64
    val::Float64
    padding::NTuple{7, Float64}
end

function check_run_time(ns_budget, start_float, end_float, t_num_floats,
                        func_name, data_format, rounding, fastmath_on,
                        search)
    # Measure the time taken to test one invocation of the function and
    # calculate the number of tests to do according to the search strategy
    # duration chosen by the user.
    start_time = time_ns()
    spawn_threads(start_float, end_float, t_num_floats,
                  func_name, data_format, rounding, fastmath_on,
                  100000, search)
    end_time = time_ns()
    elapsed_time = (end_time - start_time)/100000
    return floor(ns_budget/elapsed_time)
end

const max_error      = Ref{Vector{PaddedFloat64}}()
const max_input      = Ref{Vector{PaddedFloat64}}()
const max_output     = Ref{Vector{PaddedFloat64}}()
const max_ref_out    = Ref{Vector{PaddedFloat64}}()
const number_of_tests = Ref{Vector{PaddedFloat64}}()
const number_of_infs = Ref{Vector{PaddedFloat64}}()

function spawn_threads(start_float, end_float, t_num_floats,
                       func_name, data_format, rounding, fastmath_on,
                       tests_to_do, search)

    thread_tasks = Task[]
    for tn in 1:Threads.nthreads()

        # Calculate the sub-interval end for this thread.
        sub_int_start = Err.nextfloatn(start_float, (tn-1)*t_num_floats, data_format)
        if tn==Threads.nthreads()
            sub_int_end = end_float
        else
            sub_int_end = Err.nextfloatn(sub_int_start, t_num_floats-1, data_format)
        end

        # Run all intervals except the last in the separate threads.
        if (tn != Threads.nthreads())
            if search == "exhaustive"
                t = Threads.@spawn (max_error[][tn].val,
                                    max_input[][tn].val,
                                    max_output[][tn].val,
                                    max_ref_out[][tn].val,
                                    number_of_tests[][tn].val,
                                    number_of_infs[][tn].val) =
                                        Err.function_max_error_exhaustive(
                                            func_name, data_format, rounding, fastmath_on,
                                            sub_int_start, sub_int_end)
            else
                t = Threads.@spawn (max_error[][tn].val,
                                    max_input[][tn].val,
                                    max_output[][tn].val,
                                    max_ref_out[][tn].val,
                                    number_of_tests[][tn].val,
                                    number_of_infs[][tn].val) =
                                        Err.function_max_error_fixed_step(
                                            func_name, data_format, rounding, fastmath_on,
                                            sub_int_start, sub_int_end, tests_to_do)
            end
            push!(thread_tasks, t)
        else
            # Run the last interval in the main thread.
            if search == "exhaustive"
                (max_error[][tn].val,
                 max_input[][tn].val,
                 max_output[][tn].val,
                 max_ref_out[][tn].val,
                 number_of_tests[][tn].val,
                 number_of_infs[][tn].val) =
                     Err.function_max_error_exhaustive(
                         func_name, data_format, rounding, fastmath_on,
                         sub_int_start, sub_int_end)
            else
                (max_error[][tn].val,
                 max_input[][tn].val,
                 max_output[][tn].val,
                 max_ref_out[][tn].val,
                 number_of_tests[][tn].val,
                 number_of_infs[][tn].val) =
                     Err.function_max_error_fixed_step(
                         func_name, data_format, rounding, fastmath_on,
                         sub_int_start, sub_int_end, tests_to_do)
            end
        end
    end

    wait.(thread_tasks)

    return (max_error, max_input, max_output, max_ref_out, number_of_tests, number_of_infs)
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

    max_error[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]
    max_input[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]
    max_output[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]
    max_ref_out[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]
    number_of_tests[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]
    number_of_infs[] =
        [PaddedFloat64(0.0, ntuple(_ -> 0.0, 7)) for _ in 1: Threads.nthreads()+1]

    for task in tasks

        task_name = task.name
        data_format = task.format
        search = task.search
        rounding = task.rounding
        fastmath_on = task.fastmath
        printstyled("Running task: $task_name\n", color=:blue)

        # Set MPFR global precision to be 20 bits more than the specified format's.
        if (data_format == "binary16")
            setprecision(BigFloat, 63)
        elseif (data_format == "binary32")
            setprecision(BigFloat, 63)
        else
            setprecision(BigFloat, 127)
        end

        printstyled("Format: $data_format Search: $search Rounding: \
                    $rounding Fastmath: $fastmath_on\n", color=:green)

        # Results table formatting. Each task specied in the JSON file has
        # a result .txt file named accordingly.
        fe = FormatExpr("{1:<10s} {2:>15s} {3:>30s} {4:>30s} {5:>30s} {6:>20s} {7:>20s}\n")
        result_table_head = format(fe, "Function", "ULPs", "Input", "Output",
                                   "MPFR", "Tests", "Infs")
        fe_hex = FormatExpr("{1:<10s} {2:>15s} {3:>30s} {4:>30s} {5:>20s}\n")
        result_table_head_hex = format(fe_hex, "Function", "ULPs", "Input", "Output",
                                       "Tests")
        open(joinpath(output_dir, "$task_name.txt"), "w") do file
            write(file, result_table_head)
        end
        open(joinpath(output_dir, "HEX_$task_name.txt"), "w") do file_hex
            write(file_hex, result_table_head_hex)
        end
        fe = FormatExpr("{1:<10s} {2:>15.10f} {3:>30.15e} \
                    {4:>30.15e} {5:>30.15e} {6:>20d} {7:>20d}\n")
        fe_hex = FormatExpr("{1:<10s} {2:>15.10f} {3:>#30x} {4:>#30x} \
                    {5:>20d}\n")

        # Caluclate how many tests of this functions can be done in the
        # running time, approximately, by a single thread.
        ns_budget = 0
        if search == "seconds"
            ns_budget = 10^9
        elseif search == "minutes"
            ns_budget = 60*10^9
        elseif search == "hours"
            ns_budget = 3600*10^9
        elseif search == "days"
            ns_budget = 24*3600*10^9
        end

        # Loop through the functions list of a particular format.
        for (func_name, v) in Functions.functions_dict[data_format]

            start_float = Functions.functions_dict[data_format][func_name][1]
            end_float = Functions.functions_dict[data_format][func_name][2]

            # Check the number of floating-point numbers in the function's input
            # domain; used in calculating sub-intervals for the different threads.
            num_floats = Err.number_of_floats_in_interval(
                start_float, end_float, data_format)
            t_num_floats = floor(num_floats/Threads.nthreads())

            # The search strategy actually used for this function: a time-budgeted
            # search falls back to the exhaustive one when the budget allows testing
            # more values than there are in the domain. This must not leak into the
            # strategy of the other functions of the task.
            func_search = search
            tests_to_do = 0
            if isa(search, Integer)
                tests_to_do = search;
            elseif search != "exhaustive"
                tests_to_do = check_run_time(ns_budget, start_float, end_float, t_num_floats,
                                             func_name, data_format, rounding, fastmath_on, search)
                if tests_to_do*Threads.nthreads() > num_floats
                    func_search = "exhaustive"
                end
            end

            if func_search == "exhaustive"
                @printf("Running %d tests (search strategy exhaustive) for the function \
                         %s with %d threads \n", num_floats, func_name, Threads.nthreads())
            else
                @printf("Running %d tests (search strategy \"%s\") for the function \
                    %s with %d threads \n",
                        tests_to_do*Threads.nthreads(), func_search, func_name, Threads.nthreads())
            end
            flush(stdout)

            spawn_threads(start_float, end_float, t_num_floats,
                          func_name, data_format, rounding, fastmath_on,
                          tests_to_do, func_search)

            # Run tests on special inputs
            if (data_format != "binary16")
                input_set = Functions.spec_inputs_dict[data_format][func_name]
                (max_error[][Threads.nthreads()+1].val,
                 max_input[][Threads.nthreads()+1].val,
                 max_output[][Threads.nthreads()+1].val,
                 max_ref_out[][Threads.nthreads()+1].val,
                 number_of_tests[][Threads.nthreads()+1].val,
                 number_of_infs[][Threads.nthreads()+1].val) =
                     Err.function_max_error_special_inputs(
                         func_name, data_format, rounding, fastmath_on, input_set)
            else
                number_of_infs[][Threads.nthreads()+1].val = 0.0
            end

            # Find index of the maximum error and report to the output file.
            i = findmax([s.val for s in max_error[]])[2]
            line = format(fe, func_name, trunc(max_error[][i].val, digits=10),
                          max_input[][i].val, max_output[][i].val,
                          max_ref_out[][i].val, sum([s.val for s in number_of_tests[]]),
                          sum([s.val for s in number_of_infs[]]))
            open(joinpath(output_dir, "$task_name.txt"), "a") do file
                write(file, line);
            end
            line = format(fe_hex, func_name, trunc(max_error[][i].val, digits=10),
                          reinterpret(Err.uint_formats[data_format],
                                      Err.formats[data_format](max_input[][i].val)),
                          reinterpret(Err.uint_formats[data_format],
                                      Err.formats[data_format](max_output[][i].val)),
                          sum([s.val for s in number_of_tests[]]),
                          sum([s.val for s in number_of_infs[]]))
            open(joinpath(output_dir, "HEX_$task_name.txt"), "a") do file_hex
                write(file_hex, line);
            end
        end
    end
end
end
