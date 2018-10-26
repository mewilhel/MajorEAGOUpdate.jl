model = EAGO.Optimizer()

q = MOI.add_variables(model,3)

m = Model(with_optimizer(EAGO.Optimizer))
@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

@NLobjective(m, Min, exp(1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2)
#@NLobjective(m, Min, log(y - x ^ 2))
#@NLobjective(m, Min, log(x))

@constraint(m, x^2 + y <= 10)
@constraint(m, x + y == 10)

@NLconstraint(m, log(y - x ^ 2) <= 0)
backend = m.moi_backend
JuMP.optimize!(m)

println("Objective value: ", JuMP.objective_value(m))
println("x = ", JuMP.result_value(x))
println("y = ", JuMP.result_value(y))
