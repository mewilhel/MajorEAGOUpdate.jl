model = EAGO.Optimizer()

#model = Clp.Optimizer()

q = MOI.add_variables(model,2)

m = Model(with_optimizer(EAGO.Optimizer))
#m = Model(with_optimizer(Clp.Optimizer))
@variable(m, 1 <= x <= 3)
@variable(m, 1 <= y <= 3)

#@NLobjective(m, Min, x + y)
@objective(m, Min, x + y)
#@NLobjective(m, Min, log(y - x ^ 2))
#@NLobjective(m, Min, log(x))

@constraint(m, x + y <= 10)
@constraint(m, x - y <= 10)

#@NLconstraint(m, y >= 0)
@constraint(m, y >= 0)
JuMP.optimize!(m)
backend = m.moi_backend
debug1 = m.moi_backend.model.optimizer.Debug

#println("Objective value: ", JuMP.objective_value(m))
#println("x = ", JuMP.result_value(x))
#println("y = ", JuMP.result_value(y))
