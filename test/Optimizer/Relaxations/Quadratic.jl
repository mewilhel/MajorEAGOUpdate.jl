@testset "Quadratic Relaxations" begin
    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)

    func1 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(5.0,x[1])],
                                             [MOI.ScalarQuadraticTerm{Float64}(2.5,x[2],x[2])],2.0)
    func2 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(4.0,x[2])],
                                             [MOI.ScalarQuadraticTerm{Float64}(-2.2,x[1],x[2])],2.1)
    func3 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                             [MOI.ScalarQuadraticTerm{Float64}(2.1,x[1],x[1])],2.2)
    func4 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                             [MOI.ScalarQuadraticTerm{Float64}(-2.1,x[1],x[1])],2.2)
    func5 = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm{Float64}(3.0,x[3])],
                                             [MOI.ScalarQuadraticTerm{Float64}(2.1,x[1],x[1])],2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(model, func1, set1)
    MOI.add_constraint(model, func2, set1)
    MOI.add_constraint(model, func3, set2)
    MOI.add_constraint(model, func4, set2)
    MOI.add_constraint(model, func5, set3)

    EAGO.Quadratic_Convexity!(model)

    @test model.ConstraintConvexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(1)] == false
    @test model.ConstraintConvexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(2)]
    @test model.ConstraintConvexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(1)]
    @test model.ConstraintConvexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(2)] == false
    @test model.ConstraintConvexity[MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(1)] == false

    target = EAGO.Optimizer()
    x = MOI.add_variables(target,3)

    model.VItoSto[1] = 1; model.VItoSto[2] = 2; model.VItoSto[3] = 3;
    model.StoToVI = EAGO.ReverseDict(model.VItoSto)

    n = EAGO.NodeBB(Float64[1.0,5.0,8.0], Float64[2.0,6.0,9.0], -Inf, Inf, 2, 1, true)

    r = EAGO.DefaultRelaxationScheme()

    EAGO.RelaxQuadratic!(target, model, n, r)

    @test target.LinearLEQConstraints[1][1].terms[1].coefficient == 5.0
    @test target.LinearLEQConstraints[1][1].terms[1].variable_index.value == 1
    @test target.LinearLEQConstraints[1][1].terms[2].coefficient == 2.5
    @test target.LinearLEQConstraints[1][1].terms[2].variable_index.value == 2
    @test target.LinearLEQConstraints[1][1].constant == 1.375
    @test target.LinearLEQConstraints[1][2].upper == 1.0
    @test target.LinearLEQConstraints[1][3] == 2

    @test target.LinearLEQConstraints[2][1].terms[1].coefficient == -11.0
    @test target.LinearLEQConstraints[2][1].terms[1].variable_index.value == 1
    @test isapprox(target.LinearLEQConstraints[2][1].terms[2].coefficient, -0.4, atol=1E7)
    @test target.LinearLEQConstraints[2][1].terms[2].variable_index.value == 2
    @test target.LinearLEQConstraints[2][1].constant == 24.1
    @test target.LinearLEQConstraints[2][2].upper == 1.0
    @test target.LinearLEQConstraints[2][3] == 2

    @test isapprox(target.LinearLEQConstraints[3][1].terms[1].coefficient, -0.63, atol=1E7)
    @test target.LinearLEQConstraints[3][1].terms[1].variable_index.value == 1
    @test target.LinearLEQConstraints[3][1].terms[2].coefficient == -3.0
    @test target.LinearLEQConstraints[3][1].terms[2].variable_index.value == 3
    @test target.LinearLEQConstraints[3][1].constant == 2.0
    @test target.LinearLEQConstraints[3][2].upper == -2.0
    @test target.LinearLEQConstraints[3][3] == 2

    @test target.LinearLEQConstraints[4][1].terms[1].coefficient == 2.1
    @test target.LinearLEQConstraints[4][1].terms[1].variable_index.value == 1
    @test target.LinearLEQConstraints[4][1].terms[2].coefficient == -3.0
    @test target.LinearLEQConstraints[4][1].terms[2].variable_index.value == 3
    @test target.LinearLEQConstraints[4][1].constant == -2.725
    @test target.LinearLEQConstraints[4][2].upper == -2.0
    @test target.LinearLEQConstraints[4][3] == 2

    @test isapprox(target.LinearLEQConstraints[5][1].terms[1].coefficient, 6.3, atol=1E7)
    @test target.LinearLEQConstraints[5][1].terms[1].variable_index.value == 1
    @test target.LinearLEQConstraints[5][1].terms[2].coefficient == 3.0
    @test target.LinearLEQConstraints[5][1].terms[2].variable_index.value == 3
    @test target.LinearLEQConstraints[5][1].constant == -2.0
    @test target.LinearLEQConstraints[5][2].upper == 3.0
    @test target.LinearLEQConstraints[5][3] == 2

    @test isapprox(target.LinearLEQConstraints[6][1].terms[1].coefficient, -6.3, atol=1E7)
    @test target.LinearLEQConstraints[6][1].terms[1].variable_index.value == 1
    @test target.LinearLEQConstraints[6][1].terms[2].coefficient == -3.0
    @test target.LinearLEQConstraints[6][1].terms[2].variable_index.value == 3
    @test target.LinearLEQConstraints[6][1].constant == 2.0
    @test target.LinearLEQConstraints[6][2].upper == -3.0
    @test target.LinearLEQConstraints[6][3] == 2
end
