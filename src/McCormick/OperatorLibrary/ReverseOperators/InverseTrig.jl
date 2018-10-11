"""
Reverse `asin`.
"""
function asin_rev!(y::Interval, x::Interval)  # y = asin(x)
    h = lo(half_pi)
    y = y ∩ IntervalType(-h, h)
    x = sin(y)
end

"""
Reverse `acos`.
"""
function acos_rev!(y::MC, x::MC)
        y = y ∩ IntervalType(0.0,hi(two_pi))
        x = x ∩ cos(y)
end

"""
Reverse `atan`.
"""
function atan_rev!(y::MC, x::MC)
        y = y ∩ IntervalType(-lo(half_pi),hi(half_pi))
        x = x ∩ tan(y)
end
