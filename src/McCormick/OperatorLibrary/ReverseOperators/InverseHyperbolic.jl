"""
Reverse function for `asinh`.
"""
function asinh_rev!(y::MC,x::MC)
    x = x ∩ sinh(y)
end

"""
Reverse function for `acosh`.
"""
function acosh_rev!(y::MC,x::MC)
    y = y ∩ IntervalType(0.0,∞)
    x = x ∩ cosh(y)
end

"""
Reverse function for `atanh`.
"""
function atanh_rev!(y::MC,x::MC)
    x = x ∩ tanh(y)
end
