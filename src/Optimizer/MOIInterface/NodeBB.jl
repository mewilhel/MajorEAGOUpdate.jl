@enum Trilean TRUE FALSE EITHER

"""
    EAGO.NodeData

Stores information associated with each node in Branch & Bound tree.
"""
mutable struct NodeBB
    LowerVar::Vector{Float64}
    UpperVar::Vector{Float64}
    LowerBound::Float64
    UpperBound::Float64
    Depth::Int
    LastBranch::Int
    DirBranch::Bool
end
NodeBB() = NodeBB(Float64[],Float64[],-Inf,Inf,0,-1,false)

# Access functions for broadcasting data easily
LowerVar(x::NodeBB) = x.LowerVar
UpperVar(x::NodeBB) = x.UpperVar
LowerBound(x::NodeBB) = x.LowerBound
UpperBound(x::NodeBB) = x.UpperBound
Depth(x::NodeBB) = x.Depth
LastBranch(x::NodeBB) = x.LastBranch
DirBranch(x::NodeBB) = x.DirBranch

"""
EAGO.NodeHistory

Stores historical information associated with solving the problem.
"""
mutable struct NodeHistory
    LowerCount::Int                       # Number of Lower Bounding Problems Solved
    UpperCount::Int                       # Number of Upper Bounding Problems Solved
    LowerBound::Dict{Int,Float64}         # Interation --> Lower Bound
    UpperBound::Dict{Int,Float64}         # Interation --> Upper Bound
    LowerTime::Dict{Int,Float64}          # Interation --> Lower Problem Run Time
    UpperTime::Dict{Int,Float64}          # Interation --> Upper Problem Run Time
    PreprocessTime::Dict{Int,Float64}     # Interation --> Preprocess Problem Run Time
    PostprocessTime::Dict{Int,Float64}    # Interation --> Postprocess Problem Run Time
    CutCount::Dict{Int,Int}              # Iteration --> Number of Cuts Used
    Count::Dict{Int,Int}                  # Interation --> Number of Nodes
end
NodeHistory() = NodeHistory(0,0,Dict{Int,Float64}(0 => -Inf),
                                Dict{Int,Float64}(0 => Inf),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Int}(0 => 0),
                                Dict{Int,Int}(0 => 0))

function copy(x::NodeBB)
    return NodeBB(x.LowerVar, x.UpperVar, x.LowerBound, x.UpperBound, x.Depth, x.LastBranch, x.DirBranch)
end

function SameBox(x::NodeBB,y::NodeBB,atol::Float64)
    bool = true
    for i in 1:length(x.LowerVar)
        x.LowerVar[i] != y.LowerVar[i]

    end
    return bool
end
