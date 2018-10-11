#=
@set_cvtrait Base.abs2(x)
@set_cvtrait Base.inv(x)
=#

#                                         #Convexity       #Mono     #Switch    #Domain
@set_cvtrait Base.deg2rad(x)               :Affine      :Increasing    -Inf    -Inf  Inf
@set_cvtrait Base.rad2deg(x)               :Affine      :Increasing    -Inf    -Inf  Inf
@set_cvtrait Base.:+(x)                    :Affine      :Increasing    -Inf    -Inf  Inf
@set_cvtrait Base.transpose(x)             :Affine      :Increasing    -Inf    -Inf  Inf

@set_cvtrait Base.:-(x)                    :Affine      :Decreasing    -Inf    -Inf  Inf

@set_cvtrait Base.log(x)                   :Concave     :Increasing    -Inf    0.0   Inf
@set_cvtrait Base.log10(x)                 :Concave     :Increasing    -Inf    0.0   Inf
@set_cvtrait Base.log2(x)                  :Concave     :Increasing    -Inf    0.0   Inf
@set_cvtrait Base.log1p(x)                 :Concave     :Increasing    -Inf   -1.0   Inf
@set_cvtrait Base.sqrt(x)                  :Concave     :Increasing    -Inf    0.0   Inf
@set_cvtrait Base.acosh(x)                 :Concave     :Increasing    -Inf    1.0   Inf

@set_cvtrait Base.exp(x)                   :Convex      :Increasing    -Inf   -Inf   Inf
@set_cvtrait Base.exp2(x)                  :Convex      :Increasing    -Inf   -Inf   Inf
@set_cvtrait Base.exp10(x)                 :Convex      :Increasing    -Inf   -Inf   Inf
@set_cvtrait Base.expm1(x)                 :Convex      :Increasing    -Inf   -Inf   Inf

@set_cvtrait Base.abs(x)                   :Convex      :DecrToInrc    -Inf   -Inf   Inf
@set_cvtrait Base.cosh(x)                  :Convex      :DecrToIncr     0.0   -Inf   Inf

@set_cvtrait Base.atan(x)                  :Convexoconcave :Increasing  0.0  -Inf  Inf
@set_cvtrait Base.asinh(x)                 :Convexoconcave :Increasing  0.0  -Inf  Inf
@set_cvtrait Base.tanh(x)                  :Convexoconcave :Increasing  0.0  -Inf  Inf
@set_cvtrait SpecialFunctions.erf(x)       :Convexoconcave :Increasing  0.0  -Inf  Inf

@set_cvtrait Base.asech(x)                 :Convexoconcave :Decreasing  0.5   0.0  1.0
@set_cvtrait Base.acos(x)                  :Convexoconcave :Decreasing  0.0  -1.0  1.0
@set_cvtrait SpecialFunctions.erfcinv(x)   :Convexoconcave :Decreasing  0.0  -Inf  Inf

@set_cvtrait Base.tan(x)                   :Concavoconvex  :Increasing  0.0  -pi/2  pi/2
@set_cvtrait Base.tand(x)                  :Concavoconvex  :Increasing  0.0  -90   90
@set_cvtrait Base.cbrt(x)                  :Concavoconvex  :Increasing  0.0  -Inf  Inf
@set_cvtrait Base.sinh(x)                  :Concavoconvex  :Increasing  0.0  -Inf  Inf
@set_cvtrait Base.asin(x)                  :Concavoconvex  :Increasing  0.0  -1.0  1.0
@set_cvtrait Base.atanh(x)                 :Concavoconvex  :Increasing  0.0  -1.0  1.0
@set_cvtrait SpecialFunctions.erfinv(x)    :Concavoconvex  :Increasing  0.0  -1.0  1.0

@set_cvtrait SpecialFunctions.erfc(x)      :Concavoconvex  :Decreasing  0.0  -Inf  Inf

@set_cvtrait Base.sech(x)                  :Convexoconcave :IncrToDecr  -Inf -Inf  Inf
