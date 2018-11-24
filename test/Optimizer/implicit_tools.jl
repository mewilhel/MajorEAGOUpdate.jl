#@testset "Implicit Lower Evaluator" begin
#end

lower_eval = ImplicitLowerEvaluator{1}()
test1 = lower_eval.obj_eval == false

f_obj(x,p) = x[1]
g_constr(x,p) = [-1.0*x[1], 1.0*p[1]]
h_imp(x,p) = [x[1]^2 + p[1]*x[1] + 4]

pl = 6.0; p = 7.5; pu = 9.0; xl = -0.78; xu = -0.4;
np = 1; nx = 1; ng = 2;
user_sparse = [(1,1)]

build_lower_evaluator!(lower_eval, obj = f_obj, constr = g_constr,
                       impfun = h_imp, np = np, nx = nx, ng = ng,
                       user_sparse = user_sparse)

lower_eval.current_node = EAGO.NodeBB(Float64[pl, xl],Float64[pu, xu],-Inf,Inf,1,1,false)
lower_eval.last_node =  EAGO.NodeBB(Float64[],Float64[],-Inf,Inf,1,1,false)

EAGO.relax_implicit!(lower_eval,p)
test2 = lower_eval.ref_p
test3 = lower_eval.P
test4 = lower_eval.X
test5 = lower_eval.state_relax
test6 = lower_eval.state_ref_relaxation
test7 = lower_eval.obj_eval
test8 = lower_eval.cnstr_eval

EAGO.relax_objective!(lower_eval)
test9 = lower_eval.obj_relax
test10 = lower_eval.obj_eval == true

EAGO.relax_constraints!(lower_eval)
test11 = lower_eval.cnstr_relax
test12 = lower_eval.cnstr_eval == true

test12 = MOI.eval_objective(lower_eval, p)

g = zeros(2)
MOI.eval_constraint(lower_eval, g, p)
test13 = g

df = zeros(1)
MOI.eval_objective_gradient(lower_eval, df, p)
test14 = df

test15 = MOI.jacobian_structure(lower_eval)

dg = zeros(2,1)
MOI.eval_constraint_jacobian(lower_eval,dg,p)
test16 = dg

w = zeros(2)
y = zeros(2)
MOI.eval_constraint_jacobian_product(lower_eval, y, p, w)
test17 = y

w = zeros(2)
y = zeros(2)
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
test00 = upper_eval

y = [p; -0.6]

build_upper_evaluator!(upper_eval, obj = f_obj, constr = g_constr,
                       impfun = h_imp, np = np, nx = nx, ng = ng,
                       user_sparse = user_sparse)

EAGO.calc_functions!(upper_eval,y)
test10 = upper_eval.func_eval
test11 = upper_eval.value_storage


test12 = MOI.eval_objective(upper_eval, y)

g = zeros(4)
MOI.eval_constraint(upper_eval, g, y)
test13 = g

df = zeros(2)
MOI.eval_objective_gradient(upper_eval, df, y)
test14 = df

test15 = MOI.jacobian_structure(upper_eval)

dg = zeros(4,2)
MOI.eval_constraint_jacobian(upper_eval,dg,y)
test16 = dg

#=
w = zeros(2)
out = zeros(4)
MOI.eval_constraint_jacobian_product(upper_eval, out, y, w)
test17 = y

w = zeros(2)
out = zeros(4)
MOI.eval_constraint_jacobian_transpose_product(upper_eval, out, y, w)
test17 = y
=#
features = MOI.features_available(upper_eval)
