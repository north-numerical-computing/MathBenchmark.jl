module FileIO

using JSON

# Valid values for various fields in the JSON config file.
format_list = ["binary16", "binary32", "binary64"]
search_list = ["exhaustive", "seconds", "minutes", "hours", "days"]
rounding_list = ["RN", "RZ", "RU", "RD"]

# Read a JSON file and return a Dict of tasks.
function read_input_file(filename::String)

    # Read the json file and parse it.
    json_string = read(filename, String)
    data = JSON.parse(json_string)

    # Convert the data to a Dict.
    tasks = Dict{String, Dict{String, Any}}()
    for (task, details) in data
        d = Dict{String, Any}()
        for (k, v) in details
            d[k] = v
        end
        tasks[task] = d
    end

    return tasks
end
export read_input_file

# Validation checks for an input task.
function validate_tasks(task)

    f = 1
    s = 1
    r = 1
    fm = 0
    for (k, v) in task
        if k == "format"
            if !(v in format_list)
                printstyled("INVALID: $k\n", color=:red)
                printstyled("Invalid format: $v\n", color=:red)
                println("Valid formats: $format_list")
                return f, s, r, fm
            end
        f = v
        elseif k == "search"
            if !(v in search_list)
                printstyled("INVALID $k\n", color=:red)
                printstyled("Invalid search: $v\n", color=:red)
                println("Valid search: $search_list")
                return f, s, r, fm
            end
        s = v
        elseif k == "rounding"
            if !(v in rounding_list)
                printstyled("INVALID $k\n", color=:red)
                printstyled("Invalid rounding: $v\n", color=:red)
                println("Valid rounding: $rounding_list")
                return f, s, r, fm
            end
        r = v
        elseif k == "fastmath"
            if !(v == 1 || v == 0)
                printstyled("INVALID $k\n", color=:red)
                printstyled("Invalid fastmath flag: $v\n", color=:red)
                println("Valid fastmath: 0 or 1")
                return f, s, r, fm
            end
        fm = Bool(v)
        end
    end
    return f, s, r, fm
end


end
