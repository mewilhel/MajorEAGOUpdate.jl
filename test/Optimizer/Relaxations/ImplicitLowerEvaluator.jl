@testset "1D Implicit Lower Evaluator" begin
    low_eval =  EAGO.ImplicitLowerEvaluator{1}()
    new_node = EAGO.NodeBB(Float64[6.0,-0.78],Float64[9.0,-0.40],-Inf,Inf,0,-1,false)
    ptest = [7.5]; ptest2 = [6.5]
    EAGO.set_current_node!(low_eval, new_node)
    @test low_eval.current_node.LowerVar == Float64[6.0,-0.78]
    @test low_eval.current_node.UpperVar == Float64[9.0,-0.40]

    f(x,p) = x[1]
    g(x,p) = [x[1]]
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
    @test low_eval.objective_fun == f
    @test low_eval.constraints_fun == g
    @test low_eval.state_fun == h
    @test low_eval.np == np
    @test low_eval.nx == nx
    @test low_eval.ng == ng
    @test low_eval.imp_opts.nx == nx
    @test low_eval.imp_opts.np == np
    @test low_eval.state_relax == [zero(MC{1})]
    @test low_eval.cnstr_relax == [zero(MC{1})]
    @test low_eval.state_ref_relaxation == [[zero(MC{1})]]
    @test isnan(low_eval.last_p[1])
    @test low_eval.ref_p == [0.0]

    # initial implicit eval
    EAGO.relax_implicit!(low_eval, ptest)
    @test low_eval.obj_eval == false
    @test low_eval.cnstr_eval == false
    @test low_eval.ref_p == [7.5]
    @test isnan(low_eval.last_p[1])
    @test isapprox(low_eval.P[1].lo, 6.0, atol=1E-4) && isapprox(low_eval.P[1].hi, 9.0, atol=1E-4)
    @test isapprox(low_eval.X[1].lo, -0.78000, atol=1E-4) && isapprox(low_eval.X[1].hi, -0.4, atol=1E-4)

    test16 = low_eval.state_relax
    @test isapprox(-0.630784, low_eval.state_relax[1].cv, atol=1E-4)
    @test isapprox(-0.520912, low_eval.state_relax[1].cc, atol=1E-4)
    @test isapprox(0.0941469, low_eval.state_relax[1].cv_grad[1], atol=1E-4)
    @test isapprox(0.0455312, low_eval.state_relax[1].cc_grad[1], atol=1E-4)
    @test isapprox(-0.77205, low_eval.state_relax[1].Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_relax[1].Intv.hi, atol=1E-4)

    test17 = low_eval.state_ref_relaxation
    @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].cv, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].cc, atol=1E-4)
    @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cv_grad[1], atol=1E-4)
    @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cc_grad[1], atol=1E-4)
    @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].Intv.hi, atol=1E-4)

    # second implicit eval
    EAGO.relax_implicit!(low_eval, ptest2)
    @test low_eval.obj_eval == false
    @test low_eval.cnstr_eval == false
    @test low_eval.ref_p[1] == 7.5
    @test isnan(low_eval.last_p[1])
    @test isapprox(-0.7249310, low_eval.state_relax[1].cv, atol=1E-4)
    @test isapprox(-0.6653304, low_eval.state_relax[1].cc, atol=1E-4)
    @test isapprox(0.0941469, low_eval.state_relax[1].cv_grad[1], atol=1E-4)
    @test isapprox(0.157755, low_eval.state_relax[1].cc_grad[1], atol=1E-4)
    @test isapprox(-0.772005, low_eval.state_relax[1].Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_relax[1].Intv.hi, atol=1E-4)
    @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].cv, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].cc, atol=1E-4)
    @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cv_grad[1], atol=1E-4)
    @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cc_grad[1], atol=1E-4)
    @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].Intv.hi, atol=1E-4)

    # test objective eval
    EAGO.relax_objective!(low_eval)
    test31 = low_eval.obj_relax
    @test low_eval.obj_eval
    @test isapprox(-0.7249310, low_eval.obj_relax.cv, atol=1E-4)
    @test isapprox(-0.6653304, low_eval.obj_relax.cc, atol=1E-4)
    @test isapprox(0.0941469, low_eval.obj_relax.cv_grad[1], atol=1E-4)
    @test isapprox(0.157755, low_eval.obj_relax.cc_grad[1], atol=1E-4)
    @test isapprox(-0.772005, low_eval.obj_relax.Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.obj_relax.Intv.hi, atol=1E-4)

    # test constraint eval
    EAGO.relax_constraints!(low_eval)
    @test isapprox(-0.7249310, low_eval.cnstr_relax[1].cv, atol=1E-4)
    @test isapprox(-0.6653304, low_eval.cnstr_relax[1].cc, atol=1E-4)
    @test isapprox(0.0941469, low_eval.cnstr_relax[1].cv_grad[1], atol=1E-4)
    @test isapprox(0.157755, low_eval.cnstr_relax[1].cc_grad[1], atol=1E-4)
    @test isapprox(-0.772005, low_eval.cnstr_relax[1].Intv.lo, atol=1E-4)
    @test isapprox(-0.4, low_eval.cnstr_relax[1].Intv.hi, atol=1E-4)
    @test low_eval.cnstr_eval

    objval = MOI.eval_objective(low_eval, ptest)
    @test isapprox(-0.63078416, objval, atol=1E-4)

    gval = zeros(1)
    MOI.eval_constraint(low_eval, gval, ptest)
    @test isapprox(-0.63078416, gval[1], atol=1E-4)

    df = zeros(1)
    MOI.eval_objective_gradient(low_eval, df, ptest)
    @test isapprox(0.09414689, df[1], atol=1E-4)

    jacobian_sparsity = MOI.jacobian_structure(low_eval)
    @test jacobian_sparsity[1][1] == 1
    @test jacobian_sparsity[1][2] == 1
    @test low_eval.jacobian_sparsity[1][1] == 1
    @test low_eval.jacobian_sparsity[1][2] == 1

    # throws errors
    @test_throws ErrorException MOI.hessian_lagrangian_structure(low_eval)
    @test_throws ErrorException EAGO._hessian_lagrangian_structure(low_eval)

    # no error
    @test_nowarn MOI.initialize(low_eval, [:Grad,:Jac])

    # throws error
    @test_throws ErrorException MOI.objective_expr(low_eval)
    @test_throws ErrorException MOI.constraint_expr(low_eval)

    jac_val = zeros(1,1)
    MOI.eval_constraint_jacobian(low_eval, jac_val, ptest)
    @test isapprox(0.0941468, jac_val[1], atol=1E-4)

    w = 0.5
    jac_prod_val = zeros(1)
    MOI.eval_constraint_jacobian_product(low_eval, jac_prod_val, ptest, w)
    @test isapprox(jac_prod_val[1], 0.047073445, atol=1E-4)

    jact_prod_val = zeros(1)
    MOI.eval_constraint_jacobian_transpose_product(low_eval, jact_prod_val, ptest, w)
    @test isapprox(jact_prod_val[1], 0.047073445, atol=1E-4)

    features = MOI.features_available(low_eval)
    @test features[1] == :Grad
    @test features[2] == :Jac
end
