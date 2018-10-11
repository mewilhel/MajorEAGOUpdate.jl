import IntervalArithmetic: IntervalBox

function IntervalBox(l::Vector{T},u::Vector{T}) where {T<:AbstractFloat}
    @assert length(l) == length(u)
    len = length(l)
    IntervalBox(SVector{len,Interval{T}}([Interval(l[i],u[i]) for i=1:len])) 
end
