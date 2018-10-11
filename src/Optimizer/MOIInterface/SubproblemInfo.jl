abstract type SubProblemInfo <: Any end

mutable struct LowerInfo <: SubProblemInfo
    Feasibility::Bool
    Value::Float64
    Solution::Vector{Float64}
    LowerVarDual::Vector{Float64}
    UpperVarDual::Vector{Float64}
end
LowerInfo() = LowerInfo(true,-Inf,Float64[],Float64[],Float64[])

mutable struct UpperInfo <: SubProblemInfo
    Feasibility::Bool
    Value::Float64
    Solution::Vector{Float64}
end
UpperInfo() = UpperInfo(true,Inf,Float64[])

mutable struct PreprocessInfo <: SubProblemInfo
    Feasibility::Bool
end
PreprocessInfo() = PreprocessInfo(true)

mutable struct PostprocessInfo <: SubProblemInfo
    Feasibility::Bool
end
PostprocessInfo() = PostprocessInfo(true)
