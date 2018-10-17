@testset "Linear Relaxations" begin
    src = EAGO.Optimizer()
    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)
    x = MOI.add_variables(src,3)

    func1 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([5.0,-2.3],[x[1],x[2]]),2.0)
    func2 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([4.0,-2.2],[x[2],x[3]]),2.1)
    func3 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([3.0,-3.3],[x[1],x[3]]),2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(src, func1, set1)
    MOI.add_constraint(src, func2, set2)
    MOI.add_constraint(src, func3, set3)

    EAGO.RelaxLinear!(src,model)

    @test model.LinearLEQConstraints[1][1].constant == 2.0
    @test model.LinearGEQConstraints[1][1].constant == 2.1
    @test model.LinearEQConstraints[1][1].constant == 2.2
    @test model.LinearLEQConstraints[1][1].terms[1].coefficient == 5.0
    @test model.LinearGEQConstraints[1][1].terms[1].coefficient == 4.0
    @test model.LinearEQConstraints[1][1].terms[1].coefficient == 3.0
    @test model.LinearLEQConstraints[1][1].terms[2].coefficient == -2.3
    @test model.LinearGEQConstraints[1][1].terms[2].coefficient == -2.2
    @test model.LinearEQConstraints[1][1].terms[2].coefficient == -3.3
    @test model.LinearLEQConstraints[1][1].terms[1].variable_index.value == 1
    @test model.LinearGEQConstraints[1][1].terms[1].variable_index.value == 2
    @test model.LinearEQConstraints[1][1].terms[1].variable_index.value == 1
    @test model.LinearLEQConstraints[1][1].terms[2].variable_index.value == 2
    @test model.LinearGEQConstraints[1][1].terms[2].variable_index.value == 3
    @test model.LinearEQConstraints[1][1].terms[2].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model.LinearLEQConstraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model.LinearGEQConstraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model.LinearEQConstraints[1][2]
    @test model.LinearEQConstraints[1][3] == 2
    @test model.LinearEQConstraints[1][3] == 2
    @test model.LinearEQConstraints[1][3] == 2
end

@testset "Quadratic Relaxations" begin
end

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
