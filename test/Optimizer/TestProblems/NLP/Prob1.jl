@testset "NLP Problem #1" begin
    # Test NLP #5
    println("----- Test Example 5 -----")
    m = Model(with_optimizer(EAGO.Optimizer))

    @variable(m, 2 <= x <= 4)
    @variable(m, 1 <= y <= 2)
    @NLobjective(m, Min, 2sqrt(x) + 3sqrt(y))

    @NLconstraint(m, exp(z) >= -10)
    @NLconstraint(m, sqrt(x) <= 2)
    @NLconstraint(m, x*y >= 4)
    @NLconstraint(m, 2sqrt(x) - 3exp(z) + x*y >= 0)

    JuMP.optimize!(m)

    backend5sol = m.moi_backend.model.optimizer.variable_primal_solution
    backend5term = m.moi_backend.model.optimizer.termination_status
    backend5pstatus = m.moi_backend.model.optimizer.primal_status
    backend5dstatus = m.moi_backend.model.optimizer.dual_status
    #debug5 = m.moi_backend.model.optimizer.Debug1
end
