module Err

include("functions.jl")
using .Functions
using Printf

neg_zeros = Dict{String, Any}()
neg_zeros["binary16"] = UInt16(0x8000)
neg_zeros["binary32"] = UInt32(0x80000000)
neg_zeros["binary64"] = UInt64(0x8000000000000000)

formats = Dict{String, Any}()
formats["binary16"] = Float16
formats["binary32"] = Float32
formats["binary64"] = Float64

uint_formats = Dict{Any, Any}()
uint_formats["binary16"] = UInt16
uint_formats["binary32"] = UInt32
uint_formats["binary64"] = UInt64
uint_formats[Float16] = UInt16
uint_formats[Float32] = UInt32
uint_formats[Float64] = UInt64

subn_mask = Dict{Any, Any}()
subn_mask[Float16] = UInt16(0x03FF)
subn_mask[Float32] = UInt32(0x007FFFFF)
subn_mask[Float64] = UInt64(0x000FFFFFFFFFFFFF)

"""
Calculate the error in ulps between a floating-point number y
and a BigFloat number z, and return it as a BigFloat number.

"""
function get_ulp_error(y::Union{Float16, Float32, Float64}, z::BigFloat)

    d = abs(y - z)
    wrap = typeof(y)
    rn = convert(wrap, z)
    ulp = eps(rn)

    # Reduce ulp if a power of two was reached by rounding up.
    if abs(rn) > abs(z) && (reinterpret(Err.uint_formats[wrap], rn) &
        Err.subn_mask[wrap] == 0) && abs(rn) != floatmin(wrap)
        ulp = ulp/2
    end

    return d / ulp
end


"""
Calculate how many floating-point values are in the provided range, inclusive.

"""
function number_of_floats_in_interval(start_float, end_float, format)

    u_format = Err.uint_formats[format]

    if sign(start_float) != sign(end_float)
        return (reinterpret(u_format, abs(start_float))
                + reinterpret(u_format, abs(end_float)) + 1)
    else
        return ((max(reinterpret(u_format, start_float),
                    reinterpret(u_format, end_float)) -
                min(reinterpret(u_format, start_float),
                    reinterpret(u_format, end_float))) + 1)
    end
end


"""
Move n number of steps from the given float x.

"""
function nextfloatn(x, n, format)

    u_format = Err.uint_formats[format]
    neg_zero = Err.neg_zeros[format]

    if (signbit(x))
        x_int = reinterpret(u_format, x) - u_format(n)
    else
        x_int = reinterpret(u_format, x) + u_format(n)
    end
    y = reinterpret(typeof(x), x_int)

    if signbit(x) != signbit(y)
        if !signbit(x)
            return typeof(x)(Inf)
        else
            x_int = neg_zero - x_int
            y = reinterpret(typeof(x), x_int)
        end
    else
        return y
    end

end


"""
    resolve_function(name, fastmath_on)

Return the Julia function to test. With `fastmath_on = true` the `Base.FastMath`
variant is used when one exists (this is what `@fastmath name(x)` would call);
functions without a fast variant (`sinpi`, `cospi`, `tanpi`) are returned as is.
"""
function resolve_function(name::AbstractString, fastmath_on::Bool)
    sym = Symbol(name)
    isdefined(Base.Math, sym) || throw(ArgumentError("unknown function: $name"))
    f = getfield(Base.Math, sym)
    if fastmath_on
        fast = get(Base.FastMath.fast_op, sym, nothing)
        fast === nothing || (f = getfield(Base.FastMath, fast))
    end
    return f
end


"""
    reference_function(name)

Return the function used to compute the correctly rounded reference in BigFloat.
"""
reference_function(name::AbstractString) = getfield(Base.Math, Symbol(name))


"""
Given an input value x in one of the three floating-point formats,
calculate y, the approximation of the function `func` for that format, and
z, the high precision correctly rounded variant in BigFloat (MPFR) computed
with `func_ref`.

"""
function calculate_function(x, func, func_ref, rounding)

    # Note: here the rounding mode could be changed before calling func, but
    # Julia currently does not provide separate mathematical functions with
    # different rounding modes.
    y = func(x)
    if isinf(y)
        @warn "The Julia mathematical function has produced infinity:\
               $func. Skipping input $x."
        return (BigFloat(0.0), y, BigFloat(NaN))
    end

    bigx = BigFloat(x)
    z = func_ref(bigx)

    # Calculate the error in ulps.
    error = get_ulp_error(y, z)

    return (error, y, z)
end


"""
Go through every floating-point value in the provided range and evaluate
the maximum ulp error for the given function.

"""
function function_max_error_exhaustive(
    func, format, rounding, fastmath_on, start_float, end_float)

    max_error::BigFloat = 0.0
    max_input = 0.0;
    max_output = 0.0;
    max_ref_out::BigFloat = 0.0;
    number_of_tests = 0;
    number_of_infs = 0;

    f_format = Err.formats[format]

    x = start_float

    f = resolve_function(func, fastmath_on)
    f_ref = reference_function(func)
    while x <= end_float
        (error, y, z) = calculate_function(x, f, f_ref, rounding)
        number_of_tests = number_of_tests + 1
        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = error
            max_input = x
            max_output = y
            max_ref_out = z
        end

        x = nextfloat(x);
    end

    max_error = Float64(max_error, RoundToZero)
    return (max_error, max_input, max_output, Float64(max_ref_out), number_of_tests, number_of_infs)
end


"""
Go through floating-point values in the provided range using a fixed-sized stepping
and evaluate the maximum ulp error for the given function.

"""
function function_max_error_fixed_step(
    func, format, rounding, fastmath_on, start_float, end_float, tests_to_do)

    max_error::BigFloat = 0.0
    max_input = 0.0;
    max_output = 0.0;
    max_ref_out::BigFloat = 0.0;
    number_of_infs = 0;

    u_format = Err.uint_formats[format]
    f_format = Err.formats[format]

    x = start_float

    step_size = ceil(number_of_floats_in_interval(x, end_float, format)/tests_to_do);

    f = resolve_function(func, fastmath_on)
    f_ref = reference_function(func)
    while x <= end_float
        (error, y, z) = calculate_function(x, f, f_ref, rounding)

        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = error
            max_input = x
            max_output = y
            max_ref_out = z
        end

        x = nextfloatn(x, step_size, format)
    end

    max_error = Float64(max_error, RoundToZero)
    return (max_error, max_input, max_output, Float64(max_ref_out), tests_to_do, number_of_infs)
end


"""
Go through every floating-point value in the provided array and evaluate
the maximum ulp error for the given function.

"""
function function_max_error_special_inputs(
    func, format, rounding, fastmath_on, input_set)

    max_error::BigFloat = 0.0
    max_input = 0.0;
    max_output = 0.0;
    max_ref_out::BigFloat = 0.0;
    number_of_tests = 0;
    number_of_infs = 0;

    f_format = Err.formats[format]

    f = resolve_function(func, fastmath_on)
    f_ref = reference_function(func)
    for x in input_set
        (error, y, z) = calculate_function(x, f, f_ref, rounding)
        number_of_tests = number_of_tests + 1
        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = error
            max_input = x
            max_output = y
            max_ref_out = z
        end
    end

    max_error = Float64(max_error, RoundToZero)
    return (max_error, max_input, max_output, Float64(max_ref_out), number_of_tests, number_of_infs)
end


end
