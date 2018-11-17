#=
# Test LP #1 (Passed 1)
println("----- Test Example 1 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@NLobjective(m, Min, x + y)

@NLconstraint(m, x + y <= 10)
@NLconstraint(m, x - y <= 10)
@NLconstraint(m, y >= 0)

JuMP.optimize!(m)
#=
backend1sol = m.moi_backend.model.optimizer.variable_primal_solution
backend1term = m.moi_backend.model.optimizer.termination_status
backend1pstatus = m.moi_backend.model.optimizer.primal_status
backend1dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug1 = m.moi_backend.model.optimizer.Debug

# Test LP #2
println("----- Test Example 2 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, -3 <= x <= -1)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)

@NLobjective(m, Min, x - y + 2z)

@NLconstraint(m, x + 2y >= -10)
@NLconstraint(m, z - 2y <= 2)
@NLconstraint(m, y >= 0)

JuMP.optimize!(m)
#=
backend2sol = m.moi_backend.model.optimizer.variable_primal_solution
backend2term = m.moi_backend.model.optimizer.termination_status
backend2pstatus = m.moi_backend.model.optimizer.primal_status
backend2dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug2 = m.moi_backend.model.optimizer.Debug

# Test LP #3
println("----- Test Example 3 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, -3 <= x <= -1)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)
@variable(m, -10 <= q <= 9)

@NLobjective(m, Min, 2x - 3y + 2z)

@NLconstraint(m, x + 2y >= -10)
@NLconstraint(m, z - 2y <= 2)
@NLconstraint(m, y >= 0)
@NLconstraint(m, q-3*z-y >= 0)

JuMP.optimize!(m)
#=
backend3sol = m.moi_backend.model.optimizer.variable_primal_solution
backend3term = m.moi_backend.model.optimizer.termination_status
backend3pstatus = m.moi_backend.model.optimizer.primal_status
backend3dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug3 = m.moi_backend.model.optimizer.Debug
=#
# Test LP #4
println("----- Test Example 4 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, -3 <= x <= -1)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)
@variable(m, -10 <= q <= 9)

@NLobjective(m, Min, 2x - 3y + 2z)

@NLconstraint(m, x + 2y >= -10)
@NLconstraint(m, z - 2y <= 2)
@NLconstraint(m, y >= 4)
@NLconstraint(m, q-3*z-y >= 0)


JuMP.optimize!(m)
#=
backend4sol = m.moi_backend.model.optimizer.variable_primal_solution
backend4term = m.moi_backend.model.optimizer.termination_status
backend4pstatus = m.moi_backend.model.optimizer.primal_status
backend4dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug4a = m.moi_backend.model.optimizer.Debug1
debug4b = m.moi_backend.model.optimizer.Debug2
#=
# Test NLP #5
println("----- Test Example 5 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 2 <= x <= 7)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)

@NLobjective(m, Min, 2sqrt(x) - 3exp(z))

@NLconstraint(m, exp(z) >= -10)
@NLconstraint(m, sqrt(x) <= 2)
@NLconstraint(m, x*y >= 4)
@NLconstraint(m, 2sqrt(x) - 3exp(z) + x*y >= 0)

JuMP.optimize!(m)
#=
backend5sol = m.moi_backend.model.optimizer.variable_primal_solution
backend5term = m.moi_backend.model.optimizer.termination_status
backend5pstatus = m.moi_backend.model.optimizer.primal_status
backend5dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug5 = m.moi_backend.model.optimizer.Debug


# Test NLP #6
println("----- Test Example 6 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, -3 <= x <= -1)
@variable(m, -2 <= y <= 2)
@variable(m, 1 <= z <= 3)

@NLobjective(m, Min, 2abs(x) - 3exp(z))

@NLconstraint(m, exp(z) >= -10)
@NLconstraint(m, x^3 <= 2)
@NLconstraint(m, x*y >= 4)
@NLconstraint(m, 2x^3 - 3exp(z) + x*y >= 0)

JuMP.optimize!(m)
#=
backend6sol = m.moi_backend.model.optimizer.variable_primal_solution
backend6term = m.moi_backend.model.optimizer.termination_status
backend6pstatus = m.moi_backend.model.optimizer.primal_status
backend6dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug6 = m.moi_backend.model.optimizer.Debug
#=
# Test QP Convex #7
println("----- Test Example 7 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@objective(m, Min, x + y)

@constraint(m, x + y <= 10)
@constraint(m, x - y <= 10)
@constraint(m, y >= 0)

JuMP.optimize!(m)
#=
backend7sol = m.moi_backend.model.optimizer.variable_primal_solution
backend7term = m.moi_backend.model.optimizer.termination_status
backend7pstatus = m.moi_backend.model.optimizer.primal_status
backend7dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug7 = m.moi_backend.model.optimizer.Debug

# Test QP Convex #7 Infeasible
println("----- Test Example 8 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@objective(m, Min, x + y)

@constraint(m, x + y <= 10)
@constraint(m, x - y <= 10)
@constraint(m, y >= 0)

JuMP.optimize!(m)
#=
backend8sol = m.moi_backend.model.optimizer.variable_primal_solution
backend8term = m.moi_backend.model.optimizer.termination_status
backend8pstatus = m.moi_backend.model.optimizer.primal_status
backend8dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug8 = m.moi_backend.model.optimizer.Debug

# Test QP Nonconvex #9
println("----- Test Example 9 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@objective(m, Min, x + y)

@constraint(m, x + y <= 10)
@constraint(m, x - y <= 10)
@constraint(m, y >= 0)

JuMP.optimize!(m)
#=
backend9sol = m.moi_backend.model.optimizer.variable_primal_solution
backend9term = m.moi_backend.model.optimizer.termination_status
backend9pstatus = m.moi_backend.model.optimizer.primal_status
backend9dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug9 = m.moi_backend.model.optimizer.Debug

# Test QP Nonconvex #10 Infeasible
println("----- Test Example 10 -----")
m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))

@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@objective(m, Min, x + y)

@constraint(m, x + y <= 10)
@constraint(m, x - y <= 10)
@constraint(m, y >= 0)

JuMP.optimize!(m)
=#
#=
backend10sol = m.moi_backend.model.optimizer.variable_primal_solution
backend10term = m.moi_backend.model.optimizer.termination_status
backend10pstatus = m.moi_backend.model.optimizer.primal_status
backend10dstatus = m.moi_backend.model.optimizer.dual_status
=#
debug10 = m.moi_backend.model.optimizer.Debug
=#
