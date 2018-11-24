@testset "Test Branch Rules" begin
    B = EAGO.Optimizer()
    B.VariableNumber = 2
    B.VariableInfo = [EAGO.VariableInfo(false,1.0,false,2.0,false,false,1.5),
                      EAGO.VariableInfo(false,2.0,false,6.0,false,false,4.0)]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)
    X1,X2 = EAGO.ContinuousRelativeBisect(B,S)

    @test isapprox(X1.LowerVar[1], 1.0; atol = 1E-4)
    @test isapprox(X1.UpperVar[1], 1.5; atol = 1E-2)
    @test isapprox(X2.LowerVar[1], 1.5; atol = 1E-2)
    @test isapprox(X2.UpperVar[1], 2.0; atol = 1E-4)
#=
    X1,X2 = PseudoCostBisect(B,S)
    @test X1.LowerVars[1] == 1.5
    @test X1.UpperVars[1] ==
    @test X2.LowerVars[1] ==
    @test X2.UpperVars[1] ==

    indx,pval = PseudoCostBranch(B,S)
=#
end

@testset "Test B&B Checks" begin
    B = EAGO.Optimizer()
    B.VariableNumber = 2
    B.VariableInfo = [EAGO.VariableInfo(false,1.0,false,2.0,false,false,1.5),
                  EAGO.VariableInfo(false,2.0,false,6.0,false,false,4.0)]
    S = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)

    @test EAGO.DefaultRepeatCheck(B,S) == false

    @test EAGO.DefaultTerminationCheck(B) == false

    B.Stack[1] = S; B.IterationLimit = -1; B.CurrentIterationCount = 2
    @test EAGO.DefaultTerminationCheck(B) == false

    B.IterationLimit = 1E8; B.NodeLimit = -1;
    @test EAGO.DefaultTerminationCheck(B) == false

    B.NodeLimit = 1E8; B.CurrentLowerInfo.Value = 1.1;
    B.CurrentUpperInfo.Value = 1.1 + 1.0E-6; B.AbsoluteTolerance = 1.0E-4
    @test EAGO.DefaultTerminationCheck(B) == false

    B.AbsoluteTolerance = 1.0E-1; B.RelativeTolerance = 1.0E-12
    @test EAGO.DefaultTerminationCheck(B) == false

    B.CurrentLowerInfo.Value = -Inf;  B.CurrentUpperInfo.Value = Inf
    B.RelativeTolerance = 1.0E10;
    @test EAGO.DefaultTerminationCheck(B) == true

    B.CurrentLowerInfo.Value = 2.1;  B.CurrentUpperInfo.Value = 2.1+1E-9
    B.RelativeTolerance = 1.0E10; B.AbsoluteTolerance = 1.0E-6
    @test EAGO.DefaultConvergenceCheck(B) == true

    B.RelativeTolerance = 1.0E-6; B.AbsoluteTolerance = 1.0E10
    @test EAGO.DefaultConvergenceCheck(B) == true

    B.RelativeTolerance = 1.0E-20; B.AbsoluteTolerance = 1.0E-20
    @test EAGO.DefaultConvergenceCheck(B) == false
end

@testset "Find Lower Bound" begin
    B = EAGO.Optimizer()
    B.GlobalUpperBound = -4.5
    B.Stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.Stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.Stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    Lower = EAGO.FindLowerBound(B)

    @test Lower == -5.0
end

@testset "Test Fathom!" begin
    B = EAGO.Optimizer()
    B.GlobalUpperBound = -4.5
    B.Stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.Stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.Stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    EAGO.Fathom!(B)

    @test length(B.Stack) == 1
    @test B.Stack[2].LowerBound == -5.0
end

@testset "Node Selection" begin
    B = EAGO.Optimizer()
    B.GlobalUpperBound = -4.5
    B.Stack[1] = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    B.Stack[2] = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    B.Stack[3] = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    key,node = EAGO.NodeSelectBest!(B)

    @test key == 2
    @test node.LowerBound == -5.0
end

@testset "Node Storage" begin
    B = EAGO.Optimizer()
    y = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -4.0, 1.0, 2, 1, true)
    y1 = EAGO.NodeBB(Float64[2.0,5.0], Float64[5.0,6.0], -5.0, 4.0, 2, 1, true)
    y2 = EAGO.NodeBB(Float64[2.0,3.0], Float64[4.0,5.0], -2.0, 3.0, 2, 1, true)
    EAGO.SingleStorage!(B,y)
    @test B.MaximumNodeID == 1

    EAGO.DefaultStorage(B,y1,y2)
    @test B.Stack[1].LowerBound == -4.0
    @test B.Stack[2].LowerBound == -5.0
    @test B.Stack[3].LowerBound == -2.0
    @test B.MaximumNodeID == 3
end
