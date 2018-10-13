@testset "Node Access Functions" begin
    x = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1, true)

    @test EAGO.LowerVar(x) == Float64[1.0,5.0]
    @test EAGO.UpperVar(x) == Float64[2.0,6.0]
    @test EAGO.LowerBound(x) == -3.4
    @test EAGO.UpperBound(x) == 2.1
    @test EAGO.Depth(x) == 2
    @test EAGO.LastBranch(x) == 1
    @test EAGO.DirBranch(x) == true
end
