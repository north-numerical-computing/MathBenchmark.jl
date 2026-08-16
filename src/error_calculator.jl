module Err

const IEEEFloat = Base.IEEEFloat

# Low-level MPFR access for the allocation-free reference.
const libmpfr = Base.MPFR.libmpfr
const MPFRRoundingMode = Base.MPFR.MPFRRoundingMode
const RNDN = Base.MPFR.MPFRRoundNearest
const RNDZ = Base.MPFR.MPFRRoundToZero

# ---------------------------------------------------------------------------
# Format helpers
# ---------------------------------------------------------------------------

const FORMATS = Dict("binary16" => Float16, "binary32" => Float32, "binary64" => Float64)

"""
    float_type(format::AbstractString)

Map a format name (`"binary16"`, `"binary32"`, `"binary64"`) to the Julia type.
"""
float_type(format::AbstractString) = FORMATS[format]

"""
    reference_precision(T)

Precision (in bits) of the MPFR reference computation for format `T`. Precisions
at word boundaries (64, 128) are avoided on purpose.
"""
reference_precision(::Type{Float16}) = 63
reference_precision(::Type{Float32}) = 63
reference_precision(::Type{Float64}) = 127

# Exponent (base 2) of the smallest positive subnormal.
min_subnormal_exponent(::Type{T}) where {T<:IEEEFloat} = exponent(nextfloat(zero(T)))

# ---------------------------------------------------------------------------
# Ordinals: a bijection between the floating-point numbers of a format and a
# contiguous range of integers that preserves the ordering, so that domains,
# chunks and strides can be handled with plain integer arithmetic:
#   ordinal(-0.0) == -1, ordinal(0.0) == 0, ordinal(nextfloat(x)) == ordinal(x) + 1
# ---------------------------------------------------------------------------

"""
    ordinal(x::T) -> Int64

Position of `x` among the values of its floating-point format, such that
`ordinal(nextfloat(x)) == ordinal(x) + 1`; `-0.0` and `0.0` are distinct values
with ordinals `-1` and `0`.
"""
@inline function ordinal(x::T) where {T<:IEEEFloat}
    bits = reinterpret(Base.uinttype(T), x)
    mag = Int64(bits & ~Base.sign_mask(T))
    return signbit(x) ? -mag - 1 : mag
end

"""
    float_from_ordinal(T, k::Integer) -> T

Inverse of [`ordinal`](@ref).
"""
@inline function float_from_ordinal(::Type{T}, k::Integer) where {T<:IEEEFloat}
    U = Base.uinttype(T)
    return k >= 0 ? reinterpret(T, U(k)) : reinterpret(T, U(-k - 1) | Base.sign_mask(T))
end

"""
    number_of_floats(lo::T, hi::T) -> Int128

How many floating-point values lie in the closed interval `[lo, hi]`.
"""
number_of_floats(lo::T, hi::T) where {T<:IEEEFloat} =
    Int128(ordinal(hi)) - Int128(ordinal(lo)) + 1

# ---------------------------------------------------------------------------
# Function resolution
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# References: how the high-precision value of a function is computed.
# `new_reference(ref)` returns the (possibly stateful) object a task must use,
# `evaluate(ref, x)` returns the reference value of `x` as a BigFloat and
# `ulp_error(ref, y, z)` the error in ulps of `y` with respect to it.
# ---------------------------------------------------------------------------

"""
Reference computed by calling the BigFloat method of the function (allocating).
The precision is the current default `BigFloat` precision.
"""
struct BigFloatReference{F}
    f::F
end
new_reference(ref::BigFloatReference) = ref
@inline evaluate(ref::BigFloatReference, x::IEEEFloat) = ref.f(BigFloat(x))
@inline ulp_error(::BigFloatReference, y::IEEEFloat, z::BigFloat) = Float64(get_ulp_error(y, z), RoundToZero)

"""
Reference computed in place with the MPFR routine `mpfr_<name>` at precision
`prec` bits, without allocations: each task gets its own `MPFRWorkspace`.
"""
struct MPFRReference{name}
    prec::Int
end
function MPFRReference(name::AbstractString, prec::Integer)
    sym = Symbol(name)
    sym in MPFR_FUNCTIONS || throw(ArgumentError("no MPFR reference for function: $name"))
    return MPFRReference{sym}(prec)
end

# Julia functions with an MPFR counterpart `mpfr_<name>(r, x, rounding)`.
const MPFR_FUNCTIONS = (:acos, :acosh, :asin, :asinh, :atan, :atanh, :cbrt, :cos, :cosh,
                        :exp, :exp10, :exp2, :log, :log10, :log1p, :log2, :sin, :sinh,
                        :sqrt, :tan, :tanh, :cospi, :sinpi, :tanpi)
for f in MPFR_FUNCTIONS
    mpfr_f = QuoteNode(Symbol(:mpfr_, f))
    @eval @inline mpfr_apply!(::Val{$(QuoteNode(f))}, r::BigFloat, x::BigFloat) =
        ccall(($mpfr_f, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), r, x, RNDN)
end

struct MPFRWorkspace{name}
    x::BigFloat   # input, exact
    z::BigFloat   # reference output
    d::BigFloat   # z - y
end

new_reference(ref::MPFRReference{name}) where {name} =
    MPFRWorkspace{name}(BigFloat(precision=ref.prec), BigFloat(precision=ref.prec), BigFloat(precision=ref.prec))

@inline function evaluate(ws::MPFRWorkspace{name}, x::IEEEFloat) where {name}
    ccall((:mpfr_set_d, libmpfr), Int32, (Ref{BigFloat}, Float64, MPFRRoundingMode), ws.x, Float64(x), RNDN)
    mpfr_apply!(Val(name), ws.z, ws.x)
    return ws.z
end

# `|y - z|` is formed in MPFR (exact whenever y is within the reference precision
# of z), scaled by the power-of-two ulp of z (exact) and truncated to a Float64,
# which gives exactly `get_ulp_error(y, z)` truncated to 53 bits. Scaling before
# the conversion matters for subnormal results, where |y - z| itself may not be
# representable as a Float64.
@inline function ulp_error(ws::MPFRWorkspace, y::T, z::BigFloat) where {T<:IEEEFloat}
    ccall((:mpfr_sub_d, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Float64, MPFRRoundingMode),
          ws.d, z, Float64(y), RNDN)
    ccall((:mpfr_mul_2si, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Clong, MPFRRoundingMode),
          ws.d, ws.d, -ulp_exponent(T, z), RNDN)
    return abs(ccall((:mpfr_get_d, libmpfr), Float64, (Ref{BigFloat}, MPFRRoundingMode), ws.d, RNDZ))
end

"""
    ulp_exponent(T, z::BigFloat) -> Int

Base-2 exponent of the unit in the last place of the binade of `z` in format `T`
(the ulp is `2^ulp_exponent(T, z)`), with the spacing of the subnormals below
`floatmin(T)`. Exact powers of two get the spacing of the binade above them
(`ulp(T, 2^k) == eps(T(2^k))`), consistently with `get_ulp_error`.
"""
@inline function ulp_exponent(::Type{T}, z::BigFloat) where {T<:IEEEFloat}
    emin = min_subnormal_exponent(T)
    if ccall((:mpfr_zero_p, libmpfr), Int32, (Ref{BigFloat},), z) != 0
        return emin
    end
    # z = m * 2^e with 0.5 <= |m| < 1, so floor(log2|z|) = e - 1.
    e = Int(ccall((:mpfr_get_exp, libmpfr), Clong, (Ref{BigFloat},), z))
    return max(e - precision(T), emin)
end

"""
    ulp(T, z::BigFloat) -> Float64

Unit in the last place of the binade of `z` in format `T`, see [`ulp_exponent`](@ref).
"""
ulp(::Type{T}, z::BigFloat) where {T<:IEEEFloat} = ldexp(1.0, ulp_exponent(T, z))

"""
    get_ulp_error(y::T, z::BigFloat) -> BigFloat

Error in ulps between the computed value `y` and the high-precision reference `z`.
The ulp is that of the binade of `z` (exact powers of two get the spacing of the
binade above them, subnormals the spacing of the subnormals).
"""
function get_ulp_error(y::T, z::BigFloat) where {T<:IEEEFloat}
    d = abs(y - z)
    rn = T(z)
    u = eps(rn)
    # Reduce ulp if a power of two was reached by rounding up.
    if abs(rn) > abs(z) && (reinterpret(Base.uinttype(T), rn) & Base.significand_mask(T)) == 0 &&
       abs(rn) != floatmin(T)
        u = u / 2
    end
    return d / u
end

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

"""
Outcome of testing a set of inputs of format `T`.

* `max_error`, `input`, `output`, `reference` describe the worst case found
  (`max_error < 0` if nothing was tested);
* `ntests` is the number of inputs tested;
* `ninfs` counts inputs for which the Julia function returned an infinity;
* `nfailures` counts inputs whose result cannot be compared to the reference:
  an `Inf`/`NaN` result where the correctly rounded reference is finite (spurious
  overflow), a finite/`NaN` result where the correctly rounded reference is
  infinite (missed overflow or pole), or an infinity of the wrong sign. Failures
  do not enter `max_error`. An infinity that matches an infinite correctly
  rounded reference (overflow, pole) is only counted in `ninfs`.
"""
struct Result{T<:IEEEFloat}
    max_error::Float64
    input::T
    output::T
    reference::BigFloat
    ntests::Int
    ninfs::Int
    nfailures::Int
end

Result{T}() where {T} = Result{T}(-1.0, zero(T), zero(T), BigFloat(0.0), 0, 0, 0)

"""
    combine(a::Result, b::Result) -> Result

Merge two results: the worst case of the two (`a` on ties) and the sums of the counters.
"""
function combine(a::Result{T}, b::Result{T}) where {T}
    best = a.max_error >= b.max_error ? a : b
    return Result{T}(best.max_error, best.input, best.output, best.reference,
                     a.ntests + b.ntests, a.ninfs + b.ninfs, a.nfailures + b.nfailures)
end

# Mutable accumulator used inside the kernels; converted to a `Result` at the end.
mutable struct Accumulator{T<:IEEEFloat}
    max_error::Float64
    input::T
    output::T
    reference::BigFloat
    ntests::Int
    ninfs::Int
    nfailures::Int
end

Accumulator{T}() where {T} = Accumulator{T}(-1.0, zero(T), zero(T), BigFloat(0.0), 0, 0, 0)

Result(acc::Accumulator{T}) where {T} =
    Result{T}(acc.max_error, acc.input, acc.output, acc.reference, acc.ntests, acc.ninfs, acc.nfailures)

# Copy of a BigFloat with the same precision (BigFloat(z) would alias z).
function copy_bigfloat(z::BigFloat)
    r = BigFloat(precision=precision(z))
    ccall((:mpfr_set, libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, MPFRRoundingMode), r, z, RNDN)
    return r
end

"""
    test_input!(acc, ref, f, x)

Evaluate `f(x)` and its reference, classify the outcome (see [`Result`](@ref))
and update the accumulator `acc`.
"""
@inline function test_input!(acc::Accumulator{T}, ref, f::F, x::T) where {T<:IEEEFloat, F}
    y = f(x)
    z = evaluate(ref, x)
    zr = T(z)   # correctly rounded reference in the format
    acc.ntests += 1
    isinf(y) && (acc.ninfs += 1)
    if isfinite(y) && isfinite(zr)
        err = ulp_error(ref, y, z)
        if err > acc.max_error
            acc.max_error = err
            acc.input = x
            acc.output = y
            acc.reference = copy_bigfloat(z)
        end
    elseif !(isinf(y) && isinf(zr) && signbit(y) == signbit(zr))
        acc.nfailures += 1
    end
    return nothing
end

"""
    test_range(reference, f, ::Type{T}, k_first, k_last, step) -> Result{T}

Test `f` on the values with ordinals `k_first, k_first + step, ... <= k_last`.
"""
function test_range(reference, f::F, ::Type{T}, k_first::Integer, k_last::Integer, step::Integer) where {F, T<:IEEEFloat}
    ref = new_reference(reference)
    acc = Accumulator{T}()
    k = Int128(k_first)
    last = Int128(k_last)
    step = Int128(step)
    while k <= last
        test_input!(acc, ref, f, float_from_ordinal(T, Int64(k)))
        k += step
    end
    return Result(acc)
end

"""
    test_values(reference, f, xs) -> Result{T}

Test `f` on every value of the collection `xs`.
"""
function test_values(reference, f::F, xs::AbstractVector{T}) where {F, T<:IEEEFloat}
    ref = new_reference(reference)
    acc = Accumulator{T}()
    for x in xs
        test_input!(acc, ref, f, x)
    end
    return Result(acc)
end

end
