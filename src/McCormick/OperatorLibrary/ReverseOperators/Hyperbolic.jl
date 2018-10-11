"""
Reverse function for `sinh`.
"""
function sinh_rev!(y::MC, x::MC)
    x = x ∩ asinh(y)
end

"""
Reverse function for `cosh`.
"""
function cosh_rev!(y::MC,x::MC)
    y = y ∩ IntervalType(1.0,∞)
    x = x ∩ acosh(y)
end

"""
Reverse function for `tanh`.
"""
function tanh_rev!(y::MC,x::MC)
    y = y ∩ IntervalType(-1.0,1.0)
    x = x ∩ atanh(y)
end
