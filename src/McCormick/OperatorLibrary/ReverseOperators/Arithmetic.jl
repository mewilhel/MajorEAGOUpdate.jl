"""
Reverse minus
"""
function plus_rev!(a::MC, b::MC, c::MC)  # a = b + c
    b = b ∩ (a - c)
    c = c ∩ (a - b)
end
plus_rev!(a,b,c) = plus_rev!(promote(a,b,c)...)

"""
Reverse minus
"""
function minus_rev!(a::MC, b::MC, c::MC)  # a = b - c
    b = b ∩ (a + c)
    c = c ∩ (b - a)
end
minus_rev!(a::MC, b::MC) = (b = -a; return (a, b))     # a = -b
minus_rev!(a,b,c) = minus_rev!(promote(a,b,c)...)
minus_rev!(a,b) = minus_rev!(promote(a,b)...)

"""
Reverse multiplication
"""
function mul_rev!(a::MC, b::MC, c::MC)  # a = b * c
    ((0.0 ∉ a) ||  (0.0 ∉ b)) && (c = c ∩ (a / b))
    ((0.0 ∉ a) ||  (0.0 ∉ c)) && (b = b ∩ (a / c))
end
mul_rev!(a,b,c) = mul_rev!(promote(a,b,c)...)

"""
Reverse division
"""
function div_rev!(a::MC, b::MC, c::MC)  # a = b / c
    b = b ∩ (a * c)
    c = c ∩ (b / a)
end
div_rev!(a,b,c) = div_rev!(promote(a,b,c)...)

"""
Reverse inverse
"""
function inv_rev!(a::MC, b::MC)  # a = inv(b)
    b = b ∩ inv(a)
end
inv_rev!(a,b) = inv_rev!(promote(a,b)...)

"""
Reverse power
"""
function pow_rev!(a::MC, b::MC, c::MC)  # a = b^c
    b = b ∩ ( a^(inv(c) ))
    c = c ∩ (log(a) / log(b))
end
pow_rev!(a, b, c) = pow_rev!(promote(a, b, c)...)


"""
Reverse square root
"""
function sqrt_rev!(a::MC, b::MC)  # a = sqrt(b)
    b = b ∩ (a^2)
end
sqrt_rev(a,b) = sqrt_rev(promote(a,b)...)

"""
Reverse sqr
"""
sqr_rev(f, x)  = pow_rev!(f,x,2)


"""
Reverse abs
"""
abs_rev!(a::MC, b::MC) = ()
#=
function abs_rev!(y::MC{N}, x::MC{N}) where N   # y = abs(x); refine x

    y_new = y ∩ (0..∞)

    x1 = y_new ∩ x
    x2 = -(y_new ∩ (-x))

    cv =
    cc =
    cv_grad =
    cc_grad =
    Intv = hull(Intv(x1),Intv(x2))


    y = MC{N}(cv, cc, Intv, cv_grad, cc_grad, y.cnst)

    return
end
=#
