m = Model(with_optimizer(EAGO.Optimizer))
@variable(m, x)
@variable(m, y)

@NLobjective(m, Min, sin(1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2)

@constraint(m, x^2 + y <= 10)
@constraint(m, x + y == 10)
@constraint(m, y >= 0)

@NLconstraint(m, cos(y - x ^ 2) <= 0)
ttt = JuMP.build(m)
source_evaluator = JuMP.NLPEvaluator(m)
built_evaluator = EAGO.Build_NLP_Evaluator(MC{2},source_evaluator)

#source_evaluator = MOI.initialize(jeval, Symbol[:Grad])
#d = m.nlp_data.evaluator

#=
@testset "NLP Evaluator" begin
end
=#
