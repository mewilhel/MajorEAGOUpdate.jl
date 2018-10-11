abstract type AbstractScheme end
#=
"""
Relaxation Scheme Linear is a function that takes a linear constraint and either
adds a relaxed constraint to an NLP evaluation block or to the optimizer directly.
"""
=#
mutable struct RelaxationScheme <: AbstractScheme
    OptimizerType::Symbol
    LinearRelax::Function      # true -> NLP, false -> constraint
    LinearRelaxChk::Function
    QuadraticRelax::Function
    QuadraticRelaxChk::Function
    ExprLib::Dict{Symbol,Function}
    LinearRelaxedAfterLoad::Dict{Int,Bool}
    QuadraticRelaxedAfterLoad::Dict{Int,Bool}
    NonlinearRelaxedAfterLoad::Dict{Int,Bool}
end

macro DefineDefaultLinear(set)
    expr = quote
                DefaultLinear(f::MOI.ScalarAffineFunction{Float64},s::$set,n::NodeBB) = (f,s)
                DefaultLinearCheck(f::MOI.ScalarAffineFunction{Float64},s::$set,n::NodeBB) = false
           end
    return esc(expr)
end

@DefineDefaultLinear(MOI.GreaterThan{Float64})
@DefineDefaultLinear(MOI.LessThan{Float64})
@DefineDefaultLinear(MOI.EqualTo{Float64})

DefaultRelaxationScheme() = RelaxationScheme(:LP,DefaultLinear,DefaultLinearCheck,
                                             x->x, x->x,
                                             Dict{Symbol,Function}(),
                                             Dict{Int,Bool}(),
                                             Dict{Int,Bool}(),
                                             Dict{Int,Bool}())
