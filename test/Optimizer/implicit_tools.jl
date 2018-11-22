#@testset "Implicit Lower Evaluator" begin
#end

p = ; g = ; df = ; y = ; w = ;

lower_eval = ImplicitLowerEvaluator{3}()
test1 = lower_eval.obj_eval == false

relax_implicit!(lower_eval,p)
test2 = d.ref_p
test3 = d.P
test4 = d.X
test5 = d.state_relax
test6 = d.state_ref_relaxation
test7 = d.obj_eval
test8 = d.cnstr_eval

relax_objective!(lower_eval)
test9 = lower_eval.obj_relax
test10 = lower_eval.obj_eval == true

relax_constraints!(lower_eval)
test11 = d.cnstr_relax
test12 = d.cnstr_eval == true

test12 = MOI.eval_objective(lower_eval, p)

g = ###
MOI.eval_constraint(lower_eval, g, p)
test13 = g

df = ###
MOI.eval_objective_gradient(lower_eval, df, p)
test14 = df

test15 = MOI.jacobian_structure(lower_eval)

dg = ###
MOI.eval_constraint_jacobian(lower_eval,dg,p)
test16 = dg

w = ###
y = ###
MOI.eval_constraint_jacobian_product(lower_eval, y, p, w)
test17 = y

w = ###
y = ###
MOI.eval_constraint_jacobian_transpose_product(lower_eval, y, p, w)
test17 = y

features = MOI.features_available(lower_eval)

# No error
#=
MOI.initialize(lower_eval, Symbol[:Grad])
MOI.initialize(lower_eval, Symbol[:Hess])
=#

# Has Error
#=
MOI.objective_expr(lower_eval)
MOI.constraint_expr(lower_eval)
=#

#@testset "Implicit Upper Evaluator" begin
#end
upper_eval = ImplicitUpperEvaluator()
