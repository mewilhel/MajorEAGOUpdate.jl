#=
@testset "Begin Implicit Solver 5.1 Example" begin

    opt = EAGO.Optimizer()

    obj(x,p) = x[1]
    h(x,p) = [x[1]^2 + x[1]*p[1] + 4.0]
    hjac(x,p) = [2.0*x[1]+p[1]]

    pl = 6.0; pu = 9.0; xl = -0.78; xu = -0.4

    var, opt = SolveImplicit(f, h, xl, xu, pl, pu, opt, f = obj, hj = hjac)
    pval = JuMP.value(var[1])
    fval = JuMP.objective_value(opt)
    tstatus = JuMP.termination_status(opt)
    pstatus = JuMP.primal_status(opt)
end
=#

opt = EAGO.Optimizer()

obj(x,p) = x[1]
h(x,p) = [x[1]^2 + x[1]*p[1] + 4.0]
hjac(x,p) = [2.0*x[1]+p[1]]

pl = [6.0]; pu = [9.0]; xl = [-0.78]; xu = [-0.4]

var, opt = SolveImplicit(obj, h, xl, xu, pl, pu, opt, hj = hjac)
#=
JuMP.value(var[1])
pval = JuMP.value(var[1])
fval = JuMP.objective_value(opt)
tstatus = JuMP.termination_status(opt)
pstatus = JuMP.primal_status(opt)
=#
