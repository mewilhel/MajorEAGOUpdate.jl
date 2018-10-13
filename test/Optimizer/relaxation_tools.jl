m = Model(with_optimizer(EAGO.Optimizer))
@variable(m, x)
@variable(m, y)

@NLobjective(m, Min, sin(1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2)

@constraint(m, x^2 + y <= 10)
@constraint(m, x + y == 10)
@constraint(m, y >= 0)

@NLconstraint(m, cos(y - x ^ 2) <= 0)
#ttt = JuMP.build(m)
source_evaluator = JuMP.NLPEvaluator(m)
MOI.initialize(source_evaluator , Symbol[:Grad])
opt = m.moi_backend.model.optimizer

# Check build
built_evaluator = EAGO.Build_NLP_Evaluator(MC{2}, source_evaluator, opt)

# Add current node and define point
built_evaluator.current_node = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)
xpoint = Float64[3.0,4.0]

#=
# Tests evaluation points
EAGO.forward_eval(built_evaluator.objective.storage,
                  built_evaluator.objective.nd, built_evaluator.objective.adj,
                  built_evaluator.objective.const_values, built_evaluator.parameter_values,
                  built_evaluator.current_node, xpoint, built_evaluator.subexpression_values)
EAGO.forward_eval_all(built_evaluator, xpoint)
EAGO.reverse_eval(built_evaluator.objective.storage, built_evaluator.objective.nd, built_evaluator.objective.adj)
EAGO.reverse_eval_all(built_evaluator, xpoint)
EAGO.forward_reverse_pass(built_evaluator, xpoint)

# Test evaluation features
MOI.features_available(d::Evaluator)
MOI.eval_objective(d::Evaluator, x)
MOI.eval_constraint(d::Evaluator, g, x)
MOI.eval_objective_gradient(d::Evaluator, df, x)
MOI.jacobian_structure(d::Evaluator)
MOI.eval_constraint_jacobian(d::Evaluator)
MOI.eval_constraint_jacobian_product(d::Evaluator, y, x, w)
MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
=#

#=
@testset "NLP Evaluator" begin

#MOI.hessian_lagrangian_structure(d::Evaluator)
#MOI.eval_hessian_lagrangian_product(d::Evaluator, h, x, v, σ, μ)
#MOI.eval_hessian_lagrangian(d::Evaluator, H, x, σ, μ)
end
=#

#=
@testset "Linear Relaxations" begin
end
=#

#=
@testset "Quadratic Relaxations" begin
end
=#
