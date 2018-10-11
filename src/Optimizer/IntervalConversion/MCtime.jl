using IntervalArithmetic, BenchmarkTools, StaticArrays


import Base.asinh

a = Interval(rand(),rand()+2.0)
b = Interval(rand(),rand()+2.0)
c = Interval(-0.25,0.3)
c1 = interval(-1.2,1.3)


#const half_pi_Float64 = IntervalArithmetic.multiply_by_positive_constant(0.5, IntervalArithmetic.pi_interval(Float64))
#=
function fast_find_quadrants_tan(x::Float64)
    temp = IntervalArithmetic.atomic(Interval{Float64}, x) /half_pi_Float64

    return SVector(floor(temp.lo), floor(temp.hi))
end
=#
# Checks fast tangent
function fast_tan(a::Interval{Float64})
    isempty(a) && return a

    diam(a) > pi_interval(Float64).lo && return entireinterval(a)

    # No allocs
    lo_quadrant = minimum(fast_find_quadrants_tan(a.lo))
    hi_quadrant = maximum(fast_find_quadrants_tan(a.hi))
    # 2 allocs

    lo_quadrant_mod = mod(lo_quadrant, 2)
    hi_quadrant_mod = mod(hi_quadrant, 2)
    # 2 allocs

    if lo_quadrant_mod == 0 && hi_quadrant_mod == 1
        # check if really contains singularity:
        if hi_quadrant * half_pi_Float64 ⊆ a
            return entireinterval(a)  # crosses singularity
        end

    elseif lo_quadrant_mod == hi_quadrant_mod && hi_quadrant > lo_quadrant
        # must cross singularity
        return entireinterval(a)

    end

    # No allocs past here
    Interval((sin(interval(a.lo))/cos(interval(a.lo))).lo,
             (sin(interval(a.hi))/cos(interval(a.hi))).hi)
end

const half_intv = one(interval(0.5))
const one_intv = one(Interval{Float64})
const two_intv = one(interval(2.0))
const log2_intv = log(interval(2.0))
const log10_intv = log(interval(10.0))

fast_inv(x::Interval{Float64}) = one_intv/x
fast_exp2(x::Interval{Float64}) = exp(log2_intv*x)
fast_exp10(x::Interval{Float64}) = exp(log10_intv*x)

function fast_tanh(x::Interval)
    t1 = exp(-two_intv*x)
    t2 = exp(two_intv*x)
    pos = (one_intv-t1)/(one_intv+t1)
    neg = (t2-one_intv)/(t2+one_intv)
    return pos ∩ neg
end

# ADD Handling for wrong domain, empty domain
function fast_asinh(x::Interval{Float64})
    isempty(x) && return x
    Interval(log(interval(x.lo)+sqrt(one_intv+pow(interval(x.lo),2))),
             log(interval(x.hi)+sqrt(one_intv+pow(interval(x.hi),2))))
end

# ADD Handling for wrong domain, empty domain
function fast_acosh(x::Interval{Float64})
    isempty(x) && return x
    (x.hi <= 1.0) && return Interval(-∞, log(interval(x.hi)+sqrt(interval(x.hi)+one_intv)*sqrt(interval(x.hi)-one_intv)))
    Interval(log(interval(x.lo)+sqrt(interval(x.lo)+one_intv)*sqrt(interval(x.lo)-one_intv)),
             log(interval(x.hi)+sqrt(interval(x.hi)+one_intv)*sqrt(interval(x.hi)-one_intv)))
end

# ADD Handling for wrong domain, empty domain
function fast_atanh(x::Interval{Float64})
    isempty(x) && return x
    (x.hi >= 1.0) && return Interval(half_intv*(log(one_intv+interval(x.lo))-log(one_intv-interval(x.lo))),∞)
    (x.lo <= -1.0)  && return Interval(-∞, half_intv*(log(one_intv+interval(x.hi))-log(one_intv-interval(x.hi))))
    Interval(half_intv*(log(one_intv+interval(x.lo))-log(one_intv-interval(x.lo))),
             half_intv*(log(one_intv+interval(x.hi))-log(one_intv-interval(x.hi))))
end

asinh(x::Interval{Float64}) = fast_asinh(x)
