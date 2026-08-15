module Err

include("functions.jl")
using .Functions
using Printf

using Libdl
const libmpfr_handle = Libdl.dlopen(:libmpfr)
const dlsym_lock = ReentrantLock()

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

const MPFR_UNARY = Dict(
    "acos"  => :mpfr_acos,
    "acosh" => :mpfr_acosh,
    "asin"  => :mpfr_asin,
    "asinh" => :mpfr_asinh,
    "atan"  => :mpfr_atan,
    "atanh" => :mpfr_atanh,
    "cbrt"  => :mpfr_cbrt,
    "cos"   => :mpfr_cos,
    "cosh"  => :mpfr_cosh,
    "exp"   => :mpfr_exp,
    "exp10" => :mpfr_exp10,
    "exp2"  => :mpfr_exp2,
    "log"   => :mpfr_log,
    "log10" => :mpfr_log10,
    "log1p" => :mpfr_log1p,
    "log2"  => :mpfr_log2,
    "sin"   => :mpfr_sin,
    "sinh"  => :mpfr_sinh,
    "sqrt"  => :mpfr_sqrt,
    "tan"   => :mpfr_tan,
    "tanh"  => :mpfr_tanh,
    "cospi" => :mpfr_cospi,
    "sinpi" => :mpfr_sinpi,
    "tanpi" => :mpfr_tanpi,
)

function mpfr_fptr(name::Symbol)
    return lock(dlsym_lock) do
        Libdl.dlsym(libmpfr_handle, name)
    end
end

"""
Calculate the error in ulps between a floating-point number y
and a BigFloat number z, and return it as a BigFloat number.

"""
function get_ulp_error(y::Union{Float16, Float32, Float64}, z::BigFloat, scratch::BigFloat)

    rnd = Base.MPFR.ROUNDING_MODE[]

    wrap = typeof(y)
    rn = convert(wrap, z)
    ulp = eps(rn)

    ccall((:mpfr_abs, :libmpfr), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Int32),
          scratch, z, rnd)

    # Reduce ulp if a power of two was reached by rounding up.
    if abs(rn) > scratch && (reinterpret(Err.uint_formats[wrap], rn) &
        Err.subn_mask[wrap] == 0) && abs(rn) != floatmin(wrap)
        ulp = ulp/2
    end

    ccall((:mpfr_sub_d, :libmpfr), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Cdouble, Int32),
          scratch, z, Float64(y), rnd)
    ccall((:mpfr_abs, :libmpfr), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Int32),
          scratch, scratch, rnd)

    k = -exponent(ulp)

    ccall((:mpfr_mul_2si, :libmpfr), Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Clong, Int32),
          scratch, scratch, k, rnd)

    return scratch
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

function set!(z::BigFloat, x::AbstractFloat)
    ccall((:mpfr_set_d, :libmpfr), Int32,
          (Ref{BigFloat}, Cdouble, Int32),
          z, Float64(x), Base.MPFR.ROUNDING_MODE[])
    return z
end

function apply!(z::BigFloat, fptr::Ptr{Cvoid}, x::BigFloat)
    ccall(fptr, Int32,
          (Ref{BigFloat}, Ref{BigFloat}, Int32),
          z, x, Base.MPFR.ROUNDING_MODE[])
    return z
end

"""
Given an input value x in one of the three floating-point formats,
calculate y, the approximation of the function for that format, and
z, the high precision correctly rounded variant in BigFloat (MPFR).

"""
function calculate_function(x, func, apply_fn, rounding, fastmath_on, bigx::BigFloat, bigz::BigFloat, scratch::BigFloat)
    # Note: here the rounding mode could be changed before calling func, but
    # Julia currently does not provide separate mathematical functions with
    # different rounding modes.
    if fastmath_on
        y = @fastmath(func(x))
    else
        y = func(x)
    end

    if isinf(y)
        @warn "The Julia mathematical function has produced infinity:\
               $func. Skipping input $x."
        return (BigFloat(0.0), y, BigFloat(NaN))
    end

    set!(bigx, x)
    apply_fn(bigz, bigx)

    error = get_ulp_error(y, bigz, scratch)
    return (error, y, bigz)
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

    bigx = BigFloat(0.0)
    bigz = BigFloat(0.0)
    scratch = BigFloat(0.0)

    fptr = mpfr_fptr(MPFR_UNARY[func])
    apply_fn = (dest, x) -> apply!(dest, fptr, x)

    f = getfield(Base.Math, Symbol(func));
    while x <= end_float
        (error, y, z) = calculate_function(x, f, apply_fn, rounding,
                                           fastmath_on, bigx, bigz, scratch)

        number_of_tests = number_of_tests + 1
        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = deepcopy(error)
            max_input = x
            max_output = y
            max_ref_out = deepcopy(bigz)
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

    bigx = BigFloat(0.0)
    bigz = BigFloat(0.0)
    scratch = BigFloat(0.0)

    fptr = mpfr_fptr(MPFR_UNARY[func])
    apply_fn = (dest, x) -> apply!(dest, fptr, x)

    f = getfield(Base.Math, Symbol(func))
    while x <= end_float
        (error, y, z) = calculate_function(x, f, apply_fn, rounding,
                                           fastmath_on, bigx, bigz, scratch)

        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = deepcopy(error)
            max_input = x
            max_output = y
            max_ref_out = deepcopy(z)
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

    bigx = BigFloat(0.0)
    bigz = BigFloat(0.0)
    scratch = BigFloat(0.0)

    fptr = mpfr_fptr(MPFR_UNARY[func])
    apply_fn = (dest, x) -> apply!(dest, fptr, x)

    f = getfield(Base.Math, Symbol(func))
    for x in input_set
        (error, y, z) = calculate_function(x, f, apply_fn, rounding,
                                           fastmath_on, bigx, bigz, scratch)
        number_of_tests = number_of_tests + 1
        if (isnan(z))
            number_of_infs = number_of_infs + 1
        end

        # Update max error and corresponding values.
        if error > max_error
            max_error = deepcopy(error)
            max_input = x
            max_output = y
            max_ref_out = deepcopy(z)
        end
    end

    max_error = Float64(max_error, RoundToZero)
    return (max_error, max_input, max_output, Float64(max_ref_out), number_of_tests, number_of_infs)
end


end
