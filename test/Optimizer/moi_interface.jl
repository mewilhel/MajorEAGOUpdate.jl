@testset "Set Objective" begin
    model = EAGO.Optimizer()

    MOI.set(model, MOI.ObjectiveSense(), MOI.MinSense)
    @test model.OptimizationSense == MOI.MinSense

    MOI.set(model, MOI.ObjectiveSense(), MOI.MaxSense)
    @test model.OptimizationSense == MOI.MaxSense

    MOI.set(model, MOI.ObjectiveSense(), MOI.FeasibilitySense)
    @test model.OptimizationSense == MOI.FeasibilitySense
end

@testset "Get Termination Code " begin

    model = EAGO.Optimizer()

    # Termination Status Code Checks
    model.TerminationStatusCode = MOI.Success
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.Success

    model.TerminationStatusCode = MOI.IterationLimit
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.IterationLimit
end

@testset "Add Variable, Get Number/Index" begin

    model = EAGO.Optimizer()

    # Add variable
    MOI.add_variable(model)
    nVar = MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 1

    # Add second variable
    MOI.add_variables(model,3)
    nVar = MOI.get(model, MOI.NumberOfVariables())
    @test nVar == 4

    # Get variable indices
    indx = MOI.get(model, MOI.ListOfVariableIndices())
    @test indx ==  MOI.VariableIndex[MOI.VariableIndex(1), MOI.VariableIndex(2),
                                    MOI.VariableIndex(3), MOI.VariableIndex(4)]

    @test_nowarn EAGO.check_inbounds(model, MOI.VariableIndex(1))
    @test_throws ErrorException EAGO.check_inbounds(model,MOI.VariableIndex(6))
    #EAGO.check_inbounds((model, var::MOI.SingleVariable)
end
#=
@testset "Add Linear Constraint " begin

    model = EAGO.Optimizer()

    func1 = MOI.ScalarAffineFunction{Float64}()
    func2 = MOI.ScalarAffineFunction{Float64}()
    func3 = MOI.ScalarAffineFunction{Float64}()

    set1 = MOI.LessThan{Float64}(1.0)
    set2 = MOI.GreaterThan{Float64}(2.0)
    set3 = MOI.EqualTo{Float64}(3.0)

    MOI.add_constraint(model, func1, set1)
    MOI.add_constraint(model, func2, set2)
    MOI.add_constraint(model, func3, set3)

    term1a = model.linear_le_constraints[1][1]
    term1b = model.linear_ge_constraints[1][1]
    term1c = model.linear_eq_constraints[1][1]
    term2a = model.linear_le_constraints[1][2]
    term2b = model.linear_ge_constraints[1][2]
    term2c = model.linear_eq_constraints[1][2]
end
=#
model = EAGO.Optimizer()

func1 = MOI.ScalarAffineFunction{Float64}()
func2 = MOI.ScalarAffineFunction{Float64}()
func3 = MOI.ScalarAffineFunction{Float64}()

set1 = MOI.LessThan{Float64}(1.0)
set2 = MOI.GreaterThan{Float64}(2.0)
set3 = MOI.EqualTo{Float64}(3.0)

MOI.add_constraint(model, func1, set1)
MOI.add_constraint(model, func2, set2)
MOI.add_constraint(model, func3, set3)

term1a = model.linear_le_constraints[1][1]
term1b = model.linear_ge_constraints[1][1]
term1c = model.linear_eq_constraints[1][1]
term2a = model.linear_le_constraints[1][2]
term2b = model.linear_ge_constraints[1][2]
term2c = model.linear_eq_constraints[1][2]

#=
@testset "Empty/Isempty, EaGO Model     " begin
end

@testset "Add Quadratic Constraint      " begin
end
=#
