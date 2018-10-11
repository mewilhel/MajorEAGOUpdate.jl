"""
Reverse function for `exp`.
"""
function exp_rev!(y:::MC, x::MC)
    y = y ∩ IntervalType(0.0, ∞)
    x = x ∩ log(y)
end

"""
Reverse function for `exp2`.
"""
function exp2_rev!(y::MC, x::MC)
    y = y ∩ IntervalType(0.0, ∞)
    x = x ∩ log2(y)
end

"""
Reverse function for `exp10`.
"""
function exp10_rev!(y::MC, x::MC)
    y = y ∩ IntervalType(0.0, ∞)
    x = x ∩ log10(y)
end

"""
Reverse function for `expm1`.
"""
function expm1_rev!(y::MC, x::MC)
    y = y ∩ IntervalType(-1.0, ∞)
    x = x ∩ log1p(y)
end

"""
Reverse function for `log`: ``y = \\log(x)``
"""
function log_rev!(y::MC, x::MC)
    x = x ∩ exp(y)
end

"""
Reverse function for `log2`: ``y = \\log2(x)``
"""
function log2_rev!(y::MC, x::MC)
    x = x ∩ exp2(y)
end


"""
Reverse function for `log10`: ``y = \\log10(x)``
"""
function log10_rev!(y::MC, x::MC)
    x = x ∩ exp10(y)
end

"""
Reverse function for `log1p`: ``y = \\log1p(x)``
"""
function log1p_rev!(y::MC, x::MC)
    x = x ∩ expm1(y)
end
