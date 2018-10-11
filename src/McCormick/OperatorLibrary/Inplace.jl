# Commutative binary operators
plus!(f::MC, x::MC, y::MC) = (f = x + y; return )
mult!(f::MC, x::MC, y::MC) = (f = x * y; return )
min!(f::MC, x::MC, y::MC) = (f = min(x,y); return )
max!(f::MC, x::MC, y::MC) = (f = max(x,y); return )

for commute_op in (:plus!, :mult!, :min!, :max!)
    @eval ($commute_op)(f,x,y) = ($commute_op)(promote(f,x,y)...)
    @eval ($commute_op)(f, x) =  (f = x; return )
    @eval ($commute_op)(f, x, y, zs...) =  ($commute_op)(f, ($commute_op)(x, y, zs...))
end

# Division & subtraction
div!(f::MC, x::MC, y::MC) = (f = x / y; return )
div!(f,x,y) = div!(promote(f,x,y)...)

minus!(f::MC, x::MC, y::MC) = (f = x - y; return )
minus!(f, x, y) = minus!(promote(f,x,y)...)
minus!(f::MC, x::MC) = (f = -x; return )

# Univariate operators
InPlaceList = (:exp,:exp2,:exp10,:expm1,:log,:log2,:log10,:log1p,:sin,:cos,:tan,
               :asin,:acos,:atan,:sinh,:cosh,:tanh,:asinh,:acosh,:atanh,:abs,:sqr,:sqrt)
for i in InPlaceList
    bang_sym = Symbol(String(i)*"!")
    @eval ($bang_sym)(f::MC,x::MC) = (f = ($i)(x); return )
end

# Power operators
pow!(f,x,y) = (f = pow(x,y); return )
