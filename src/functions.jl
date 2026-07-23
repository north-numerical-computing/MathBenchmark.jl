module Functions

binary16_functions = Dict{String, Tuple{Float16, Float16}}()
binary32_functions = Dict{String, Tuple{Float32, Float32}}()
binary64_functions = Dict{String, Tuple{Float64, Float64}}()

binary32_functions_spec_inputs = Dict{String, Vector{Float32}}()
binary64_functions_spec_inputs = Dict{String, Vector{Float64}}()

max_binary16 = Float16(65504.0);
max_binary32 = prevfloat(typemax(Float32));
max_binary64 = prevfloat(typemax(Float64));

functions_dict = Dict{String, Dict{}}()
functions_dict["binary16"] = binary16_functions
functions_dict["binary32"] = binary32_functions
functions_dict["binary64"] = binary64_functions

spec_inputs_dict = Dict{String, Dict{}}()
spec_inputs_dict["binary32"] = binary32_functions_spec_inputs
spec_inputs_dict["binary64"] = binary64_functions_spec_inputs

# Exact input domain intervals for 24 univariate functions
# available in the Base Julia.

# Note: functions that are not tested yet and are not in the Julia base
#    erf, erfc, expm1, j0, j1, lgamma, tgamma, y0,
#    y1, acospi, asinpi, atanpi, rsqrt

binary16_functions["acos"]  = (-1.0, 1)
binary16_functions["acosh"] = (1.0, max_binary16)
binary16_functions["asin"]  = (-1.0, 1.0)
binary16_functions["asinh"] = (-max_binary16, max_binary16)
binary16_functions["atan"]  = (-max_binary16, max_binary16)
binary16_functions["atanh"] = (nextfloat(Float16(-1.0)), prevfloat(Float16(1.0)))
binary16_functions["cbrt"]  = (-max_binary16, max_binary16)
binary16_functions["cos"]   = (-max_binary16, max_binary16)
binary16_functions["cosh"]  = (-11.78125, 11.78125)
binary16_functions["exp"]   = (-16.625, 11.0859375)
binary16_functions["exp10"] = (-7.22265625, 4.8125)
binary16_functions["exp2"]  = (-24, 15.9921875)
binary16_functions["log"]   = (2^(-24), max_binary16)
binary16_functions["log10"] = (2^(-24), max_binary16)
binary16_functions["log1p"] = (nextfloat(Float16(-1)), max_binary16)
binary16_functions["log2"]  = (2^(-24), max_binary16)
binary16_functions["sin"]   = (-max_binary16, max_binary16)
binary16_functions["sinh"]  = (-11.78125, 11.78125)
binary16_functions["sqrt"]  = (0.0, max_binary16)
binary16_functions["tan"]   = (-max_binary16, max_binary16)
binary16_functions["tanh"]  = (-4.50390625, 4.50390625)
binary16_functions["cospi"] = (-max_binary16, max_binary16)
binary16_functions["sinpi"] = (-max_binary16, max_binary16)
binary16_functions["tanpi"] = (-max_binary16, max_binary16)

binary32_functions["acos"]  = (-1, 1)
binary32_functions["acosh"] = (1.0, max_binary32)
binary32_functions["asin"]  = (-1.0, 1.0)
binary32_functions["asinh"] = (-max_binary32, max_binary32)
binary32_functions["atan"]  = (-max_binary32, max_binary32)
binary32_functions["atanh"] = (nextfloat(Float32(-1.0)), prevfloat(Float32(1.0)))
binary32_functions["cbrt"]  = (-max_binary32, max_binary32)
binary32_functions["cos"]   = (-max_binary32, max_binary32)
binary32_functions["cosh"]  = (-89.415985107421875, 89.415985107421875)
binary32_functions["exp"]   = (-103.27892303466796875, 88.72283172607421875)
binary32_functions["exp10"] = (-44.853466033935546875, 38.531841278076171875)
binary32_functions["exp2"]  = (-149, 127.99999237060546875)
binary32_functions["log"]   = (2^(-149), max_binary32)
binary32_functions["log10"] = (2^(-149), max_binary32)
binary32_functions["log1p"] = (nextfloat(Float32(-1)), max_binary32)
binary32_functions["log2"]  = (2^(-149), max_binary32)
binary32_functions["sin"]   = (-max_binary32, max_binary32)
binary32_functions["sinh"]  = (-89.415985107421875, 89.415985107421875)
binary32_functions["sqrt"]  = (0.0, max_binary32)
binary32_functions["tan"]   = (-max_binary32, max_binary32)
binary32_functions["tanh"]  = (-9.01091289520263671875, 9.01091289520263671875)
binary32_functions["cospi"] = (-max_binary32, max_binary32)
binary32_functions["sinpi"] = (-max_binary32, max_binary32)
binary32_functions["tanpi"] = (-max_binary32, max_binary32)

# Inputs causing maximum errors in AMD LibM 5.1 and Newlib 4.5.0 as reported
# in https://gitlab.inria.fr/zimmerma/math_accuracy/-/blob/master/latex/glibc242.tex?ref_type=heads
binary32_functions_spec_inputs["acos"]  = parse.(Float32, ["-0x1.737196p-2", "-0x1.0a77f4p-1"])
binary32_functions_spec_inputs["acosh"] = parse.(Float32, ["0x1.0098d6p+13", "0x1.01cb98p+0"])
binary32_functions_spec_inputs["asin"]  = parse.(Float32, ["0x1.001328p-1", "0x1.048506p-1"])
binary32_functions_spec_inputs["asinh"] = parse.(Float32, ["0x1.fffbdep+0", "0x1.0949p-1"])
binary32_functions_spec_inputs["atan"]  = parse.(Float32, ["0x1.bfcf24p-2", "0x1.626772p-1"])
binary32_functions_spec_inputs["atanh"] = parse.(Float32, ["0x1.d921aep-2", "0x1.ec2dd8p-4"])
binary32_functions_spec_inputs["cbrt"]  = parse.(Float32, ["0x1.01ff68p-127", "0x1.c828b4p-127"])
binary32_functions_spec_inputs["cos"]   = parse.(Float32, ["0x1.6d336ep-1", "0x1.921fcp+0"])
binary32_functions_spec_inputs["cosh"]  = parse.(Float32, ["0x1.4f58dap+0", "0x1.65a55ap+6"])
binary32_functions_spec_inputs["exp"]   = parse.(Float32, ["-0x1.638ffap-8", "-0x1.2e35b6p+5"])
binary32_functions_spec_inputs["exp10"] = parse.(Float32, ["-0x1.35a566p-9", "-0x1.66d3e8p+5"])
binary32_functions_spec_inputs["exp2"]  = parse.(Float32, ["-0x1.fe49bep-8", "-0x1.2a0002p+7"])
binary32_functions_spec_inputs["log"]   = parse.(Float32, ["0x1.105602p+0", "0x1.472fdep+0"])
binary32_functions_spec_inputs["log10"] = parse.(Float32, ["0x1.109116p+0", "0x1.52846ep+0"])
binary32_functions_spec_inputs["log1p"] = parse.(Float32, ["0x1.6351d8p+95", "0x1.a827b4p-2"])
binary32_functions_spec_inputs["log2"]  = parse.(Float32, ["0x1.15fde4p+0", "0x1.69cf3ep+0"])
binary32_functions_spec_inputs["sin"]   = parse.(Float32, ["0x1.48776cp+35", "0x1.ffe2a4p+101"])
binary32_functions_spec_inputs["sinh"]  = parse.(Float32, ["0x1.65a55ap+6"])
binary32_functions_spec_inputs["sqrt"]  = parse.(Float32, [])
binary32_functions_spec_inputs["tan"]   = parse.(Float32, ["0x1.ff953ep-8", "0x1.921fcp+0"])
binary32_functions_spec_inputs["tanh"]  = parse.(Float32, ["0x1.01b232p+0", "0x1.ddca18p-3"])
binary32_functions_spec_inputs["cospi"] = parse.(Float32, ["0x1.d73b56p-2"])
binary32_functions_spec_inputs["sinpi"] = parse.(Float32, ["0x1.45fcc2p-7"])
binary32_functions_spec_inputs["tanpi"] = parse.(Float32, ["0x1.fe513ep-3"])

binary64_functions["acos"]  = (-1.0, 1.0)
binary64_functions["acosh"] = (1.0, max_binary64)
binary64_functions["asin"]  = (-1.0, 1.0)
binary64_functions["asinh"] = (-max_binary64, max_binary64)
binary64_functions["atan"]  = (-max_binary64, max_binary64)
binary64_functions["atanh"] = (nextfloat(Float64(-1.0)), prevfloat(1.0))
binary64_functions["cbrt"]  = (-max_binary64, max_binary64)
binary64_functions["cos"]   = (-max_binary64, max_binary64)
binary64_functions["cosh"]  = (-710.4758600739438634263933636248111724853515625, 710.4758600739438634263933636248111724853515625)
binary64_functions["exp"]   = (-744.4400719213812180896638892590999603271484375, 709.7827128933839730962063185870647430419921875)
binary64_functions["exp10"] = (-323.30621534311575260289828293025493621826171875, 308.2547155599166899264673702418804168701171875)
binary64_functions["exp2"]  = (-1074, 1023.9999999999998863131622783839702606201171875)
binary64_functions["log"]   = (2^(-1074), max_binary64)
binary64_functions["log10"] = (2^(-1074), max_binary64)
binary64_functions["log1p"] = (nextfloat(Float64(-1)), max_binary64)
binary64_functions["log2"]  = (2^(-1074), max_binary64)
binary64_functions["sin"]   = (-max_binary64, max_binary64)
binary64_functions["sinh"]  = (-710.4758600739438634263933636248111724853515625, 710.4758600739438634263933636248111724853515625)
binary64_functions["sqrt"]  = (0.0, max_binary64)
binary64_functions["tan"]   = (-max_binary64, max_binary64)
binary64_functions["tanh"]  = (-19.06154746539849753617090755142271518707275390625, 19.06154746539849753617090755142271518707275390625)
binary64_functions["cospi"] = (-max_binary64, max_binary64)
binary64_functions["sinpi"] = (-max_binary64, max_binary64)
binary64_functions["tanpi"] = (-max_binary64, max_binary64)

# Inputs causing maximum errors in GNU libc 2.41 and IML 2025.0.0 as reported
# in https://inria.hal.science/hal-03141101v8/document
# and some arbitrarily selected worst cases from
# https://gitlab.inria.fr/core-math/core-math/-/raw/master/src/binary64/
# (.wc files for each function).
binary64_functions_spec_inputs["acos"]  = parse.(Float64, ["0x1.dffffb3488a4p-1", "0x1.6c05eb219ec46p-1",
                                                           "0x1.ffffe7507938p-1", "0x1.ffffe78d96c87p-1"])
binary64_functions_spec_inputs["acosh"] = parse.(Float64, ["0x1.0001ff6afc4bap+0", "0x1.01825ca7da7e5p+0",
                                                           "0x1.92acc901d133fp+3", "0x1.93e3865814b01p+3"])
binary64_functions_spec_inputs["asin"]  = parse.(Float64, ["-0x1.0000045b2c904p-3", "0x1.6c042a6378102p-1",
                                                           "0x1.a6a58d55e3071p-25", "0x1.a6a58d55e3072p-25"])
binary64_functions_spec_inputs["asinh"] = parse.(Float64, ["-0x1.02657ff36d5f3p-2", "-0x1.000276d9cf31ap-4"])
binary64_functions_spec_inputs["atan"]  = parse.(Float64, ["0x1.f9004c4fef9eap-4", "-0x1.ffff8020d3d1dp-7",
                                                           "0x1.f9a746fadc626p+1", "0x1.fd5ca4d9c850ep+1"])
binary64_functions_spec_inputs["atanh"] = parse.(Float64, ["-0x1.ebb5133a9d9a4p-4", "-0x1.e2cfb2667f17ep-9",
                                                           "0x1.d12ed0af1a2a4p-27", "0x1.d12ed0af1a2a5p-27"])
binary64_functions_spec_inputs["cbrt"]  = parse.(Float64, ["0x1.7a337e1ba1ec2p-257", "-0x1.f7af4893d1d51p-616",
                                                           "0x1.0b0b055fab24cp-1", "0x1.0b6b16786e828p-1"])
binary64_functions_spec_inputs["cos"]   = parse.(Float64, ["-0x1.7120161c92674p+0", "-0x1.d19ebc5567dcdp+311",
                                                           "0x1.054f8fb7ecd2bp-10", "0x1.05871c6ef5f9ap-4"])
binary64_functions_spec_inputs["cosh"]  = parse.(Float64, ["-0x1.633c654fee2bap+9", "-0x1.5a364e6b98134p+9",
                                                           "0x1.200698cf881a2p-1", "0x1.200f9eeddc5a7p-11"])
binary64_functions_spec_inputs["exp"]   = parse.(Float64, ["-0x1.49f33ad2c1c58p+9", "0x1.fce66609f7428p+5",
                                                           "0x1.0000000000068p-52", "0x1.0000000000068p-53"])
binary64_functions_spec_inputs["exp10"] = parse.(Float64, ["-0x1.57449153f316ep-7", "-0x1.5cd9d94d49a85p+1",
                                                           "0x1.14d54004e8937p-2", "0x1.14ef2dec8c0d8p-18"])
binary64_functions_spec_inputs["exp2"]  = parse.(Float64, ["-0x1.1a4ce073ea908p-5", "-0x1.8002f29666b99p-6",
                                                           "0x1.423b81681e9bap+9", "0x1.428p+9"])
binary64_functions_spec_inputs["log"]   = parse.(Float64, ["0x1.1211bef8f68e9p+0", "0x1.008000db2e8bep+0",
                                                           "0x1.dcf0ee0466b06p+5", "0x1.0d7ad3a32b788p+6"])
binary64_functions_spec_inputs["log10"] = parse.(Float64, ["0x1.de02157073b31p-1", "0x1.feda7b62c1033p-1",
                                                           "0x1.5563f85c8b3bep-997", "0x1.5853b748bc157p-997"])
binary64_functions_spec_inputs["log1p"] = parse.(Float64, ["-0x1.2bf183e0344b2p-2", "0x1.000aee2a2757fp-9",
                                                           "0x1.1f8aada39e276p-6", "0x1.1f93142cb3154p-5"])
binary64_functions_spec_inputs["log2"]  = parse.(Float64, ["0x1.1406d79e1b574p+0", "0x1.b4ebe40c95a01p+0",
                                                           "0x1.fbe3592b80acdp+0", "0x1.fc187e351a8cap+0"])
binary64_functions_spec_inputs["sin"]   = parse.(Float64, ["-0x1.f8b791cafcdefp+4", "-0x1.0e16eb809a35dp+944",
                                                           "0x1.2bf085ea1cc22p-18", "0x1.2bf60953039e1p-4"])
binary64_functions_spec_inputs["sinh"]  = parse.(Float64, ["-0x1.633c654fee2bap+9", "-0x1.adc135eb544c1p-2",
                                                           "0x1.11cee83f97d7dp+4", "0x1.120276d078a45p-17"])
binary64_functions_spec_inputs["sqrt"]  = parse.(Float64, ["0x1.fffffffffffffp-1", "0x1.fffffffffffffp-1"])
binary64_functions_spec_inputs["tan"]   = parse.(Float64, ["-0x1.317cd745dd37cp+9", "0x1.49adfd996a81dp+18",
                                                           "0x1.0484a223bae5ap-18", "0x1.049587001c4a9p-23"])
binary64_functions_spec_inputs["tanh"]  = parse.(Float64, ["0x1.e126eee514cbcp-3", "0x1.002629fd74484p+0",
                                                           "0x1.47b9905a89da1p-8", "0x1.482d9b57e7a8cp-15"])
binary64_functions_spec_inputs["cospi"] = parse.(Float64, ["-0x1.1a0a2fa299b92p+6", "-0x1.f01d619f61e5bp-8",
                                                           "0x1.508934ccfabe7p-5", "0x1.4e4d895ac7df2p-5"])
binary64_functions_spec_inputs["sinpi"] = parse.(Float64, ["-0x1.45f3e53e1d707p-7", "0x0.07f16ec91c164p-1022",
                                                           "0x1.391b1c1cfaae5p-49", "0x1.39979a8cf1447p-49"])
binary64_functions_spec_inputs["tanpi"] = parse.(Float64, ["-0x1.fae7d0ef22d4ep-2", "0x1.49b79692667bp+46",
                                                           "0x1.d5a05afa3bfcap-54", "0x1.d54d5c04ecedep-54"])


end
