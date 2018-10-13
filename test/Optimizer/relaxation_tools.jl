
src = EAGO.Optimizer()
trg = Ipopt.Optimizer()

MOI.add_variables(src,3)

@constraint(m, x + y == 10)
@constraint(m, x + 2*z <= 10)
@constraint(m, 3*y + z >= 10)

src = m.moi_backend.model.optimizer

EAGO.RelaxLinear!(src,trg)

term1a = trg.linear_le_constraints[1][1]
term1b = trg.linear_ge_constraints[1][1]
term1c = trg.linear_eq_constraints[1][1]

term2a = trg.linear_le_constraints[1][2]
term2b = trg.linear_ge_constraints[1][2]
term2c = trg.linear_eq_constraints[1][2]

#=
@testset "Linear Relaxations" begin
    linear_le_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}}
    linear_ge_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}}
    linear_eq_constraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}}
end
=#

#=
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

user_operators = built_evaluator.m.nlp_data.user_operators
nlp_data = built_evaluator.m.nlp_data
user_input_buffer = built_evaluator.jac_storage

EAGO.forward_eval(built_evaluator.objective.storage,
                  built_evaluator.objective.nd, built_evaluator.objective.adj,
                  built_evaluator.objective.const_values, built_evaluator.parameter_values,
                  built_evaluator.current_node, xpoint, built_evaluator.subexpression_values,
                  user_input_buffer, user_operators = user_operators)
=#
#EAGO.forward_eval_all(built_evaluator, xpoint)
#=
# Tests evaluation points
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
@testset "Quadratic Relaxations" begin
end
=#
