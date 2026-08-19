module Config

using JSON

# Valid values for the fields of a task in the JSON config file.
const FORMATS = ("binary16", "binary32", "binary64")
# Alternative names of the formats, mapped to the canonical name.
const FORMAT_ALIASES = Dict(
    "Float16" => "binary16",
    "half"   => "binary16",
    "Float32" => "binary32",
    "single" => "binary32",
    "Float64" => "binary64",
    "double" => "binary64",
)
const SEARCHES = ("exhaustive", "seconds", "minutes", "hours", "days")
const ROUNDINGS = ("RN", "RZ", "RU", "RD")

"""
A validated testing task from the JSON config file.

* `format`: `"binary16"`, `"binary32"` or `"binary64"` (the aliases `"Float16"`,
  `"half"`, `"Float32"`, `"single"`, `"Float64"` and `"double"` are accepted in
  the config file and mapped to these names);
* `search`: `"exhaustive"`, a time budget (`"seconds"`, `"minutes"`, `"hours"`,
  `"days"`) or a positive integer number of tests per thread;
* `rounding`: `"RN"`, `"RZ"`, `"RU"` or `"RD"` (currently without effect, Julia
  does not provide mathematical functions with different rounding modes);
* `fastmath`: whether to test the `@fastmath` variants of the functions.
"""
struct BenchmarkTask
    name::String
    format::String
    search::Union{String, Int}
    rounding::String
    fastmath::Bool
end

_invalid(name, msg) = throw(ArgumentError("task \"$name\": $msg"))

function _field(name, details, key)
    haskey(details, key) || _invalid(name, "missing field \"$key\"")
    return details[key]
end

"""
    validate_task(name, details) -> BenchmarkTask

Check the fields of a task read from the config file and throw an `ArgumentError`
with an explanatory message if any of them is missing or invalid.
"""
function validate_task(name::AbstractString, details)
    details isa AbstractDict ||
        _invalid(name, "expected an object with the fields \"format\", \"search\", \"rounding\" and \"fastmath\"")

    format = _field(name, details, "format")
    format isa AbstractString && (format = get(FORMAT_ALIASES, format, format))
    format in FORMATS ||
        _invalid(name, "invalid format $(repr(format)), valid formats: $(join(FORMATS, ", ")) \
                        or the aliases $(join(sort!(collect(keys(FORMAT_ALIASES))), ", "))")

    search = _field(name, details, "search")
    (search in SEARCHES || (search isa Integer && !(search isa Bool) && search > 0)) ||
        _invalid(name, "invalid search $(repr(search)), valid searches: $(join(SEARCHES, ", ")) or a positive integer")

    rounding = _field(name, details, "rounding")
    rounding in ROUNDINGS ||
        _invalid(name, "invalid rounding $(repr(rounding)), valid roundings: $(join(ROUNDINGS, ", "))")

    fastmath = _field(name, details, "fastmath")
    (fastmath isa Bool || (fastmath isa Integer && fastmath in (0, 1))) ||
        _invalid(name, "invalid fastmath flag $(repr(fastmath)), valid values: 0 or 1")

    return BenchmarkTask(String(name), String(format),
                         search isa Integer ? Int(search) : String(search),
                         String(rounding), Bool(fastmath))
end

"""
    read_config(filename) -> Vector{BenchmarkTask}

Read and validate the tasks of a JSON config file, in the order they appear in it.
"""
function read_config(filename::AbstractString)
    data = JSON.parsefile(filename)
    data isa AbstractDict ||
        throw(ArgumentError("$filename: expected a JSON object mapping task names to task descriptions"))
    return [validate_task(String(name), details) for (name, details) in data]
end

end
