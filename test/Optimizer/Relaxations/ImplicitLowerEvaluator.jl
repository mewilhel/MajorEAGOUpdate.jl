low_eval =  EAGO.ImplicitLowerEvaluator{1}()
new_node = EAGO.NodeBB(Float64[6.0,-0.78],Float64[9.0,-0.40],-Inf,Inf,0,-1,false)
ptest = [7.5]; ptest2 = [6.5]
EAGO.set_current_node!(low_eval, new_node)
test0A_1a = low_eval.current_node.LowerVar == Float64[6.0,-0.78]
test0A_1b = low_eval.current_node.UpperVar == Float64[9.0,-0.40]

f(x,p) = x[1]
g(x,p) = x[1]
function h(x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t4 + t3
    return [t5]
end

hj(x,p) = [2*x[1]+p[1]]
np = 1; nx = 1; ng = 1
sparse_pattern = Tuple{Int64,Int64}[(1,1)]

EAGO.build_lower_evaluator!(low_eval, h, np, nx)
EAGO.build_lower_evaluator!(low_eval, h, np, nx, obj = f)
EAGO.build_lower_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
EAGO.build_lower_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
EAGO.build_lower_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj)
EAGO.build_lower_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj, user_sparse = sparse_pattern)
#test01_1 = low_eval.objective_fun == f
#test01_2 = low_eval.constraints_fun == g
#test01_3 = low_eval.state_fun == h
#test01_4 = low_eval.np == np
#test01_5 = low_eval.nx == nx
#test01_6 = low_eval.ng == ng
#test01_7 = low_eval.imp_opts.nx == nx
#test01_8 = low_eval.imp_opts.np == np
#test01_9 = low_eval.state_relax == [zero(MC{1})]
#test01_10 = low_eval.cnstr_relax == [zero(MC{1})]
#test01_11 = low_eval.state_ref_relaxation == [[zero(MC{1})],[zero(MC{1})],[zero(MC{1})]]
#test01_14 = isnan(low_eval.last_p[1])
#test01_15 = low_eval.ref_p == [0.0]

# initial implicit eval
println("Evaluate Implicit Internal #1")
EAGO.relax_implicit!(low_eval, ptest)
test11 = low_eval.obj_eval == false
test12 = low_eval.cnstr_eval == false
test13 = low_eval.ref_p == [7.5]
test13b = low_eval.last_p == [7.5]
test14 = isapprox(low_eval.P[1].lo, 6.0, atol=1E-4) && isapprox(low_eval.P[1].hi, 9.0, atol=1E-4)
test15 = isapprox(low_eval.X[1].lo, -0.78000, atol=1E-4) && isapprox(low_eval.X[1].hi, -0.4, atol=1E-4)
test16 = low_eval.state_relax
test17 = low_eval.state_ref_relaxation

# second implicit eval
println("Evaluate Implicit Internal #2")
EAGO.relax_implicit!(low_eval, ptest2)
test11a = low_eval.obj_eval == false
test12a = low_eval.cnstr_eval == false
test13a = low_eval.ref_p == [7.5]
test13a = low_eval.last_p == [7.5]
test14a = low_eval.state_relax
#test15 = isapprox(low_eval.X[1].lo, -0.78000, atol=1E-4) && isapprox(low_eval.X[1].hi, -0.4, atol=1E-4)

#=
# test objective eval
EAGO.relax_objective!(low_eval)
test31 = low_eval.obj_relax
test32 = low_eval.obj_eval

# test constraint eval
EAGO.relax_constraints!(low_eval)
test41 = low_eval.cnstr_relax
test42 = low_eval.cnstr_eval

objval = MOI.eval_objective(low_eval, ptest)
test51 = objval

gval = zeros(1)
MOI.eval_constraint(low_eval, gval, ptest)
test61 = gval

df = zeros(1)
MOI.eval_objective_gradient(low_eval, df, ptest)
test71 = df

test81 = low_eval.jacobian_sparsity
jacobian_sparsity = MOI.jacobian_structure(low_eval)
test82 = jacobian_sparsity
test83 = low_eval.jacobian_sparsity

# throws errors
# MOI.hessian_lagrangian_structure(low_eval)
# EAGO._hessian_lagrangian_structure(low_eval)

jac_val = zeros(1)
MOI.eval_constraint_jacobian(low_eval, jac_val, ptest)
test91 = jac_val

w = 0.5
jac_prod_val = zeros(1,1)
MOI.eval_constraint_jacobian_product(low_eval, jac_prod_val, ptest, w)
test10_1 = jac_prod_val

jact_prod_val = zeros(1,1)
MOI.eval_constraint_jacobian_transpose_product(low_eval, jact_prod_val, ptest, w)
test11_1 = jact_prod_val

#
features = MOI.features_available(low_eval)
test12_1 = features

# no error
#MOI.initialize(low_eval, [:Grad,:Jac])

# throws error
# MOI.objective_expr(low_eval)
# MOI.constraint_expr(low_eval)
=#
