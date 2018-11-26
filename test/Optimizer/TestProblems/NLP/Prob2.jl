@testset "NLP Problem #2" begin
    # Test NLP #6
    println("----- Test Example 6 -----")
    m = Model(with_optimizer(EAGO.Optimizer))

    @variable(m, -3 <= x <= -1)
    @variable(m, -2 <= y <= 2)
    @variable(m, 1 <= z <= 3)

    @NLobjective(m, Min, 2abs(x) - 3exp(z))

    @NLconstraint(m, exp(z) >= -10)
    @NLconstraint(m, x^3 <= 2)
    @NLconstraint(m, x*y >= 4)
    @NLconstraint(m, 2x^3 - 3exp(z) + x*y >= 0)

    JuMP.optimize!(m)

    backend6sol = m.moi_backend.model.optimizer.variable_primal_solution
    backend6term = m.moi_backend.model.optimizer.termination_status
    backend6pstatus = m.moi_backend.model.optimizer.primal_status
    backend6dstatus = m.moi_backend.model.optimizer.dual_status

    debug6 = m.moi_backend.model.optimizer.Debug1
end
