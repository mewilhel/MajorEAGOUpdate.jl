divide_positive_constant(α, x::Interval) = @round(α*x.hi, α*x.lo)

pi_to_degree(::Type{Float64}) = divide_positive_constant(180.0, pi_interval(Float64))
pi_to_degree(::Type{T}) where {T} = 180.0/pi_interval(T)
pi_to_degree(x::T) where {T<:AbstractFloat} = PiToDeg(T)

degree_to_pi(::Type{Float64}) = multiply_by_positive_constant(1.0/180.0, pi_interval(Float64))
degree_to_pi(::Type{T}) where {T} = pi_interval(T)/180.0
degree_to_pi(x::T) where {T<:AbstractFloat} = DegToPi(T)

sind(x::Interval{T}) where T = sin(degree_to_pi(T)*x)
cosd(x::Interval{T}) where T = cos(degree_to_pi(T)*x)
tand(x::Interval{T}) where T = tan(degree_to_pi(T)*x)
secd(x::Interval{T}) where T = sec(degree_to_pi(T)*x)
cscd(x::Interval{T}) where T = csc(degree_to_pi(T)*x)
cotd(x::Interval{T}) where T = cot(degree_to_pi(T)*x)

asind(x::Interval{T}) where T = asin(pi_to_degree(T)*x)
acosd(x::Interval{T}) where T = acos(pi_to_degree(T)*x)
atand(x::Interval{T}) where T = atan(pi_to_degree(T)*x)
asecd(x::Interval{T}) where T = asec(pi_to_degree(T)*x)
acscd(x::Interval{T}) where T = acsc(pi_to_degree(T)*x)
acotd(x::Interval{T}) where T = acot(pi_to_degree(T)*x)

sinpi(x::Interval{T}) where T = sin(pi_interval(T)*x)
cospi(x::Interval{T}) where T = cos(pi_interval(T)*x)

deg2rad(x::Interval{T}) where T = degree_to_pi(T)*x
rad2deg(x::Interval{T}) where T = pi_to_degree(T)*x

transpose(x::Interval) = x
