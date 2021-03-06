"Storage for EAGO's NLP evaluator"
mutable struct FunctionSetStorage{T}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{T}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    dependent_subexpressions::Vector{Int}
end

FunctionSetStorage(T) = FunctionSetStorage{T}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],T[],Float64[],
                                           Bool[],Int[],Int[],Int[],Int[])

mutable struct SubexpressionSetStorage{T}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{T}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    linearity::JuMP.Derivatives.Linearity
end

function SubexpressionSetStorage(T,nd::Vector{JuMP.NodeData}, const_values, num_variables, subexpression_linearity, moi_index_to_consecutive_index)

    nd = JuMP.replace_moi_variables(nd, moi_index_to_consecutive_index)
    len_nd = length(nd)
    adj = adjmat(nd)
    setstorage = zeros(T,len_nd)
    numberstorage = zeros(len_nd)
    numvalued = zeros(Bool,len_nd)
    linearity = JuMP.classify_linearity(nd, adj, subexpression_linearity)

    return SubexpressionSetStorage{T}(nd, adj, const_values, setstorage, numberstorage,
                                      numvalued, linearity[1])
end

"Container for calculating relaxations of nonlinear terms"
mutable struct Evaluator{T<:Real} <: MOI.AbstractNLPEvaluator
    m::Model
    has_user_mv_operator::Bool
    parameter_values::Vector{Float64}
    variable_number::Int
    current_node::NodeBB
    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    has_reverse::Bool
    fw_repeats::Int
    fw_atol::Float64
    objective::FunctionSetStorage
    constraints::Vector{FunctionSetStorage{T}}
    subexpressions::Vector{SubexpressionSetStorage{T}}
    subexpression_order::Vector{Int}
    subexpression_values_flt::Vector{Float64}
    subexpression_values_set::Vector{T}
    subexpression_linearity::Vector{JuMP.Derivatives.Linearity}
    subexpressions_as_julia_expressions::Vector{Any}
    last_x::Vector{Float64}
    last_obj::T
    jac_storage::Vector{T}
    user_output_buffer::Vector{T}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
    function Evaluator{T}(m) where T<:Real
        d = new()
        d.m = m
        d.constraints = FunctionSetStorage[]
        d.objective = FunctionSetStorage(T)
        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0
        return d
    end
end

function set_current_node!(x::Evaluator,n::NodeBB)
    x.current_node = n
end

eltype(x::Evaluator{T}) where T  = T

include("Univariate.jl")
include("Passes.jl")
include("GetInfo.jl")
include("Load.jl")
