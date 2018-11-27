jumpmodel9 = Model(with_optimizer(EAGO.Optimizer, InitialRelaxedOptimizer = CPLEX.Optimizer()))
@variable(jumpmodel9, -5 <= x9 <= 5)
@variable(jumpmodel9, -5 <= y9 <= 5)
@NLobjective(jumpmodel9, Min, 2*x9^2-1.05*x9^4+(x9^6)/6+x9*y9+y9^2)
status9 = JuMP.optimize!(jumpmodel9)

#=
@testset "NLP Problem #3" begin
    jumpmodel9 = Model(with_optimizer(EAGO.Optimizer))
    @variable(jumpmodel9, -5 <= x9 <= 5)
    @variable(jumpmodel9, -5 <= y9 <= 5)
    @NLobjective(jumpmodel9, Min, 2*x9^2-1.05*x9^4+(x9^6)/6+x9*y9+y9^2)
    status9 = JuMP.optimize!(jumpmodel9)
    @test isapprox(getvalue(x9),0.0,atol=1E0)
    @test isapprox(getvalue(y9),0.0,atol=1E0)
    @test isapprox(getobjectivevalue(jumpmodel9),0.0,atol=1E-1)
    @test status9 == :Optimal
end
=#
