@testset "NLP Problem #8: ex4_1_3 (global library)" begin

    m = Model(with_optimizer(EAGO.Optimizer))


    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1]
    @variable(m, x[x_Idx])
    setlowerbound(x[1], 0.0)
    setupperbound(x[1], 10.0)

    # ----- Constraints ----- #
    @NLconstraint(m, e1, -(8.9248e-5*x[1]-0.0218343* (x[1])^2+0.998266* (x[1])^3-1.6995* (x[1])^4+0.2* (x[1])^5)+objvar == 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(fval,-443.6717,atol=1E-3)
    @test status_term == MOI.Success
    @test status_prim == MOI.FeasiblePoint
end
