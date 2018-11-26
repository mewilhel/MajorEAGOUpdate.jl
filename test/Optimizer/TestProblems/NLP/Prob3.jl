@testset "NLP Problem #3" begin
    jumpmodel4 = Model(with_optimizer(EAGO.Optimizer))
    @variable(jumpmodel4, -200 <= x <= -100)
    @variable(jumpmodel4, 200 <= y <= 400)
    @constraint(jumpmodel4, -500 <= x+2y <= 400)
    @NLobjective(jumpmodel4, Min, x*y)
    status4 = solve(jumpmodel4)

    @test status4 == :Optimal
    @test isapprox(getvalue(x),-200.0,atol=1E-6)
    @test isapprox(getvalue(y),300.0,atol=1E-6)
    @test isapprox(getobjectivevalue(jumpmodel4),-60000.00119999499,atol=2.0)
end
