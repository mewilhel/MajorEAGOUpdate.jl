@testset "Linear Relaxations" begin
    src = EAGO.Optimizer()
    model = EAGO.Optimizer()

    x = MOI.add_variables(model,3)
    x = MOI.add_variables(src,3)

    func1 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([5.0,-2.3],[x[1],x[2]]),2.0)
    func2 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([4.0,-2.2],[x[2],x[3]]),2.1)
    func3 = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm.([3.0,-3.3],[x[1],x[3]]),2.2)

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(src, func1, set1)
    MOI.add_constraint(src, func2, set2)
    MOI.add_constraint(src, func3, set3)

    EAGO.RelaxLinear!(src,model)

    @test model.LinearLEQConstraints[1][1].constant == 2.0
    @test model.LinearGEQConstraints[1][1].constant == 2.1
    @test model.LinearEQConstraints[1][1].constant == 2.2
    @test model.LinearLEQConstraints[1][1].terms[1].coefficient == 5.0
    @test model.LinearGEQConstraints[1][1].terms[1].coefficient == 4.0
    @test model.LinearEQConstraints[1][1].terms[1].coefficient == 3.0
    @test model.LinearLEQConstraints[1][1].terms[2].coefficient == -2.3
    @test model.LinearGEQConstraints[1][1].terms[2].coefficient == -2.2
    @test model.LinearEQConstraints[1][1].terms[2].coefficient == -3.3
    @test model.LinearLEQConstraints[1][1].terms[1].variable_index.value == 1
    @test model.LinearGEQConstraints[1][1].terms[1].variable_index.value == 2
    @test model.LinearEQConstraints[1][1].terms[1].variable_index.value == 1
    @test model.LinearLEQConstraints[1][1].terms[2].variable_index.value == 2
    @test model.LinearGEQConstraints[1][1].terms[2].variable_index.value == 3
    @test model.LinearEQConstraints[1][1].terms[2].variable_index.value == 3
    @test MOI.LessThan{Float64}(1.0) == model.LinearLEQConstraints[1][2]
    @test MOI.GreaterThan{Float64}(2.0) == model.LinearGEQConstraints[1][2]
    @test MOI.EqualTo{Float64}(3.0) == model.LinearEQConstraints[1][2]
    @test model.LinearEQConstraints[1][3] == 2
    @test model.LinearEQConstraints[1][3] == 2
    @test model.LinearEQConstraints[1][3] == 2
end
