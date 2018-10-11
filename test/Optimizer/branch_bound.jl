#=
@testset "Test Branch Rules" begin
    nbb = NodeBB()
    opt = EAGO.Optimizer()

    X1,X2 = ContinuousRelativeBisect(B,S)
    @test X1.LowerVars[2] ==
    @test X1.UpperVars[2] ==
    @test X2.LowerVars[2] ==
    @test X2.UpperVars[2] ==

    X1,X2 = PseudoCostBisect(B,S)
    @test X1.LowerVars[1] ==
    @test X1.UpperVars[1] ==
    @test X2.LowerVars[1] ==
    @test X2.UpperVars[1] ==

    indx,pval = PseudoCostBranch(B,S)
end

@testset "Test B&B Checks" begin
end

@testset "Test Fathom!" begin
end

@testset "Find Lower Bound" begin
end

@testset "Node Selection" begin
end

@testset "Node Storage" begin
end
=#
