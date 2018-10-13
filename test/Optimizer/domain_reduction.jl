@testset "Duality-Based Bound Tightening" begin

    ylower = [1.0, 1.0, 1.0, 1.0]
    yupper = [4.0, 4.0, 4.0, 4.0]
    ymult_lo = [50, 1.0, 2.0, 3.0]
    ymult_hi = [0, 1.0, 2.0, 3.0]
    yLBD = 1.0
    yUBD = 3.0
    n = EAGO.NodeBB(ylower,yupper,yLBD,yUBD,3,2,false)

    EAGO.Variable_DBBT!(n,ymult_lo,ymult_hi,yLBD,yUBD)
    @test 3.95999-1E-4 <= n.LowerVar[1] <= 3.95999+1E-4
    @test 4.0-1E-4 <= n.UpperVar[1] <= 4.0+1E-4
    @test 2.0-1E-4 <= n.LowerVar[2] <= 2.0+1E-4
    @test 4.0-1E-4 <= n.UpperVar[2] <= 4.0+1E-4
    @test 3.0-1E-4 <= n.LowerVar[3] <= 3.0+1E-4
    @test 4.0-1E-4 <= n.UpperVar[3] <= 4.0+1E-4
    @test 3.33333-1E-4 <= n.LowerVar[4] <= 3.33333+1E-4
    @test 4.0-1E-4 <= n.UpperVar[4] <= 4.0+1E-4
end

@testset "Poor Man's LP" begin
    # Puranik 2017 example
    opt1 = EAGO.Optimizer()
    x = MOI.add_variables(opt1,3)

    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[1]), MOI.GreaterThan(-2.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.LessThan(4.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[2]), MOI.GreaterThan(0.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.LessThan(1.0))
    MOI.add_constraint(opt1,MOI.SingleVariable(x[3]), MOI.GreaterThan(-1.0))

    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], -2.0), MOI.GreaterThan(2.0))
    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[2]), MOI.ScalarAffineTerm(1.0, x[3])], -1.0), MOI.LessThan(0.0))

    n1 = EAGO.NodeBB()
    n1.LowerVar = [-2.0, 0.0, -1.0]
    n1.UpperVar = [4.0, 4.0, 1.0]
    n1.LowerBound = -Inf
    n1.UpperBound = Inf
    n1.Depth = 3

    opt1.VItoSto = Dict{Int,Int}(1 => 1, 2 => 2, 3 => 3)
    opt1.NonlinearVariable = Dict{Int,Int}(1 => true, 2 => true, 3 => true)

    feas1 = EAGO.PoorManLP(opt1,n1)
    feas1 = EAGO.PoorManLP(opt1,n1)
    @test (n1.LowerVar == [2.0, 0.0, -1.0])
    @test (n1.UpperVar == [4.0, 2.0, 1.0])
    @test (feas1 == true)
end

#=
@testset "Quadratic Classification" begin
    m = EAGO.Optimizer()
end

@testset "Quadratic Domain Reduction (Univariate)" begin
    m = EAGO.Optimizer()

    n2 = NodeData()
    n2.LowerVar = [-10.0, -10.0, -10.0]
    n2.UpperVar = [10.0, 10.0, 10.0]
    n2.LowerBound = -Inf
    n2.UpperBound = Inf
    n2.Depth = 1
    a1, b1, c1 = 2.0, –4.0, 3.0
    a2, b2, c2 = 2.0, –4.0, 3.0
    a3, b3, c3 = 2.0, –4.0, 3.0
    push!(m.UniQuadraticGEQConstraints,(a1,b1,c1,1))
    push!(m.UniQuadraticLEQConstraints,(a2,b2,c2,2))
    push!(m.UniQuadraticEQConstraints,(a3,b3,c3,3))
    feas = UnivariateQuadratic(m,n2)

    #@test feas ==
    #@test n2.LowerVar[1] ==
    #@test n2.LowerVar[2] ==
    #@test n2.LowerVar[3] ==
    #@test n2.UpperVar[1] ==
    #@test n2.UpperVar[2] ==
    #@test n2.UpperVar[3] ==
end

@testset "Interval CSP" begin
    # Vigerske 2017 example
    opt1 = EAGO.Optimizer()
    x = MOI.add_variable(opt1)
    y = MOI.add_variable(opt1)
    z = MOI.add_variable(opt1)

    MOI.add_constraint(opt1,MOI.SingleVariable(x), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(x), MOI.GreaterThan(-Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(y), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(y), MOI.GreaterThan(-Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(z), MOI.LessThan(Inf))
    MOI.add_constraint(opt1,MOI.SingleVariable(z), MOI.GreaterThan(-Inf))

    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], 0.0), MOI.GreaterThan(4.0))
    MOI.add_constraint(opt1,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]), MOI.ScalarAffineTerm(1.0, x[2])], 0.0), MOI.LessThan(1.0))

    n1 = NodeData()
    n1.LowerVar = [-Inf, -Inf, 0.0]
    n1.UpperVar = [Inf, Inf, 1.0]
    n1.LowerBound = -Inf
    n1.UpperBound = Inf
    n1.Depth = 3

    # C1 = exp(x^2 + y^2) - z <= 0
    # C2 = z*sqrt(x^2+y^2) <= 1

    feas1 = EAGO.IntervalCSP(opt1,n1)
    @test (n1.LowerVar[1] == 0.0)
    @test (n1.LowerVar[2] == 0.0)
    @test (n1.LowerVar[3] == 1.0)
    @test (n1.UpperVar[1] == 0.0)
    @test (n1.UpperVar[2] == 0.0)
    @test (n1.UpperVar[3] == 1.0)
    @test feas1 == true
end

@testset "Optimization-Based Bound Tightening" begin
    # Puranik 2017 example
    opt2 = EAGO.Optimizer()
    y = MOI.add_variables(opt2,2)

    MOI.add_constraint(opt2,MOI.SingleVariable(y[1]), MOI.LessThan(5.0))
    MOI.add_constraint(opt2,MOI.SingleVariable(y[1]), MOI.GreaterThan(-3.0))
    MOI.add_constraint(opt2,MOI.SingleVariable(y[2]), MOI.LessThan(5.0))
    MOI.add_constraint(opt2,MOI.SingleVariable(y[2]), MOI.GreaterThan(-3.0))

    MOI.add_constraint(opt2,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, y[1]), MOI.ScalarAffineTerm(1.0, y[2])], 0.0), MOI.GreaterThan(0.0))
    MOI.add_constraint(opt2,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, y[1]), MOI.ScalarAffineTerm(1.0, y[2])], 0.0), MOI.LessThan(4.0))
    MOI.add_constraint(opt2,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(-1.0, y[1]), MOI.ScalarAffineTerm(1.0, y[2])], 0.0), MOI.GreaterThan(-2.0))
    MOI.add_constraint(opt2,MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(-1.0, y[1]), MOI.ScalarAffineTerm(1.0, y[2])], 0.0), MOI.LessThan(2.0))

    n2 = NodeData()
    n2.LowerVar = [-3.0, -3.0]
    n2.UpperVar = [5.0, 5.0]
    n2.LowerBound = -Inf
    n2.UpperBound = Inf
    n2.Depth = 3

    feas2 =  EAGO.OBBT(opt2,n2)
    @test (n2.LowerVar == [-1.0, -1.0])
    @test (n2.UpperVar == [3.0, 3.0])
    @test feas2 == true
end
=#
