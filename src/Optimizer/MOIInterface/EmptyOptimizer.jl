#=
Defines empty structures
=#

mutable struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
end
EmptyNLPEvaluator() = EmptyNLPEvaluator(NodeBB())

MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end


empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)


mutable struct DummyOptimizer <: MOI.AbstractOptimizer end
