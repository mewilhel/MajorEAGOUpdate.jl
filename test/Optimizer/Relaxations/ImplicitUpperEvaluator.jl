@testset "Implicit Upper Evaluator" begin

    upper_eval =  EAGO.ImplicitUpperEvaluator()
    new_node = EAGO.NodeBB(Float64[5.0, 5.0, 0.0, 0.0],
                           Float64[10.0, 10.0, 10.0, 10.0],
                           -Inf,Inf,0,-1,false)
    ptest = [7.5]; ptest2 = [6.5]
    EAGO.set_current_node!(upper_eval, new_node)

    f(x,p) = x[1]^2 + x[2]^2 + p[1]^2
    g(x,p) = [x[1]
              x[2]]
    h(x,p) = [x[1]^2 + p[1]*x[1] + 4.0
              x[2]^2 + p[2]*x[2] + 4.0]
    np = 2; nx = 2; ng = 2
    y = [6.0, 7.1, 3.2, 4.3]

    out1 = zeros(7); EAGO.reform_1!(out1,y,f,g,h,np,ng,nx)
    out2 = zeros(6); EAGO.reform_2!(out2,y,g,h,np,ng,nx)
    out3 = zeros(5); EAGO.reform_3!(out3,y,f,h,np,nx)
    out4 = zeros(4); EAGO.reform_4!(out4,y,h,np,nx)

    @test isapprox(out1[1], 64.73, atol=1E-3)
    @test isapprox(out1[2], 3.20, atol=1E-3)
    @test isapprox(out1[3], 4.30, atol=1E-3)
    @test isapprox(out1[4], 33.44, atol=1E-3)
    @test isapprox(out1[5], 53.01999, atol=1E-3)
    @test isapprox(out1[6], -33.4400, atol=1E-3)
    @test isapprox(out1[7], -53.01999, atol=1E-3)

    @test isapprox(out2[1], 3.20, atol=1E-3)
    @test isapprox(out2[2], 4.30, atol=1E-3)
    @test isapprox(out2[3], 33.44, atol=1E-3)
    @test isapprox(out2[4], 53.01999, atol=1E-3)
    @test isapprox(out2[5], -33.44, atol=1E-3)
    @test isapprox(out2[6], -53.01999, atol=1E-3)

    @test isapprox(out3[1], 64.73, atol=1E-3)
    @test isapprox(out3[2], 33.44, atol=1E-3)
    @test isapprox(out3[3], 53.01999, atol=1E-3)
    @test isapprox(out3[4], -33.44, atol=1E-3)
    @test isapprox(out3[5], -53.01999, atol=1E-3)

    @test isapprox(out4[1], 33.44, atol=1E-3)
    @test isapprox(out4[2], 53.01999, atol=1E-3)
    @test isapprox(out4[3], -33.44, atol=1E-3)
    @test isapprox(out4[4], -53.01999, atol=1E-3)

    EAGO.build_upper_evaluator!(upper_eval, h, np, nx)
    @test upper_eval.nx == 2
    @test upper_eval.ny == 4
    @test upper_eval.ng == 0
    @test upper_eval.np == 2
    @test isapprox(upper_eval.value_storage[1],0.0,atol=1E-3)
    @test isapprox(upper_eval.value_storage[2],0.0,atol=1E-3)
    @test isapprox(upper_eval.value_storage[3],0.0,atol=1E-3)
    @test isapprox(upper_eval.value_storage[4],0.0,atol=1E-3)
    @test upper_eval.diff_result == fill(0.0,(4,4))
    @test upper_eval.last_y[1] == 0.0
    @test upper_eval.last_y[2] == 0.0
    @test upper_eval.last_y[3] == 0.0
    @test upper_eval.last_y[4] == 0.0
    @test upper_eval.has_nlobj == false

    EAGO.calc_functions!(upper_eval, y)
    test_9d = upper_eval.diff_result
    test_10d = upper_eval.value_storage

    EAGO.build_upper_evaluator!(upper_eval, h, np, nx; obj = f)
    @test upper_eval.nx == 2
    @test upper_eval.ny == 4
    @test upper_eval.ng == 0
    @test upper_eval.np == 2
    test_5e = upper_eval.value_storage
    test_6e = upper_eval.diff_result
    test_7e = upper_eval.last_y
    @test upper_eval.has_nlobj

    EAGO.calc_functions!(upper_eval, y)
    test_9e = upper_eval.diff_result
    test_10e = upper_eval.value_storage
    test_11e = upper_eval.last_y

    EAGO.build_upper_evaluator!(upper_eval, h, np, nx; obj = f, constr = g, ng = 2)
    @test upper_eval.nx == 2
    @test upper_eval.ny == 4
    @test upper_eval.ng == 2
    @test upper_eval.np == 2
    test_5f = upper_eval.value_storage
    test_6f = upper_eval.diff_result
    test_7f = upper_eval.last_y
    @test upper_eval.has_nlobj

    EAGO.calc_functions!(upper_eval, y)
    test_9f = upper_eval.diff_result
    test_10f = upper_eval.value_storage

    fval = MOI.eval_objective(upper_eval, y)
    @test isapprox(fval, 64.73, atol=1E-3)

    gval = zeros(6)
    MOI.eval_constraint(upper_eval, gval, y)
    @test isapprox(gval[1], 3.20, atol=1E-3)
    @test isapprox(gval[2], 4.30, atol=1E-3)
    @test isapprox(gval[3], 33.44, atol=1E-3)
    @test isapprox(gval[4], 53.01999, atol=1E-3)
    @test isapprox(gval[5], -33.44, atol=1E-3)
    @test isapprox(gval[6], -53.01999, atol=1E-3)

    dfval = zeros(4)
    MOI.eval_objective_gradient(upper_eval, dfval, y)
    @test isapprox(dfval[1], 12.0, atol=1E-3)
    @test isapprox(dfval[2], 0.0, atol=1E-3)
    @test isapprox(dfval[3], 6.4, atol=1E-3)
    @test isapprox(dfval[4], 8.6, atol=1E-3)

    jac_struct = MOI.jacobian_structure(upper_eval)
    @test jac_struct[1][1] == 1
    @test jac_struct[1][2] == 1
    @test jac_struct[2][1] == 1
    @test jac_struct[2][2] == 2
    @test jac_struct[3][1] == 2
    @test jac_struct[3][2] == 1
    @test jac_struct[4][1] == 2
    @test jac_struct[4][2] == 2

    feat_sym = MOI.features_available(upper_eval)
    @test feat_sym[1] == :Grad
    @test feat_sym[2] == :Jac

    dg = zeros(6,4)
    MOI.eval_constraint_jacobian(upper_eval, dg, y)
    @test dg ==  [0.0  0.0   1.0   0.0;
                      0.0  0.0   0.0   1.0;
                      3.2  0.0  12.4   0.0;
                      0.0  4.3   0.0  15.7;
                     -3.2  0.0 -12.4   0.0;
                      0.0 -4.3   0.0 -15.7]

    #=
    out = zeros(6); w = fill(0.5,(6,))
    MOI.eval_constraint_jacobian_product(upper_eval, out, y, w)
    MOI.eval_constraint_jacobian_transpose_product(upper_eval, y, p, w)
    =#

    # should error
    # MOI.hessian_lagrangian_structure(upper_eval)
    # _hessian_lagrangian_structure(upper_eval)
    # MOI.objective_expr(upper_eval)
    # MOI.constraint_expr(upper_eval)
end
