@testset "Test Exponentials" begin

   mctol = 1E-4
   m = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0),
             seedg(Float64,1,2), seedg(Float64,1,2), false)

   y = McCormick.exp(m)
   yref = MC{2}(15.154262241479262, 37.30486063158251, EAGO.IntervalType(2.71828, 54.5982), SVector{2,Float64}([0.0, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.exp2(m)
   yref = MC{2}(4.0, 11.33333333333333, EAGO.IntervalType(1.99999, 16.0001), SVector{2,Float64}([2.77259, 0.0]), SVector{2,Float64}([4.66667, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.exp10(m)
   yref = MC{2}(1000.0, 6670.0000000000055, EAGO.IntervalType(9.999999999999999999, 10000.00000000001), SVector{2,Float64}([2302.5850929940457, 0.0]), SVector{2,Float64}([3330.0, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = McCormick.expm1(m)
   yref = MC{2}(6.38905609893065, 36.304860631582514, EAGO.IntervalType(1.71828, 53.5982), SVector{2,Float64}([7.38906, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

end

m = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0),
          seedg(Float64,1,2), seedg(Float64,1,2), false)

ylog = log(m)
ylog2 = log2(m)
ylog10 = log10(m)
ylog1p = log1p(m)

@testset "Test Logarithms" begin

   mctol = 1E-4
   m = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0),
             seedg(Float64,1,2), seedg(Float64,1,2), false)

   y = log(m)
   yref = MC{2}(0.9241962407465939, 0.6931471805599453, EAGO.IntervalType(0, 1.3863), SVector{2,Float64}([0.462098, 0.0]), SVector{2,Float64}([0.5, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log2(m)
   yref = MC{2}(1.3333333333333333, 1.0, EAGO.IntervalType(0, 2), SVector{2,Float64}([0.666667, 0.0]), SVector{2,Float64}([0.721348, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log10(m)
   yref = MC{2}(0.4013733275519749, 0.3010299956639812, EAGO.IntervalType(0, 0.60206), SVector{2,Float64}([0.200687, 0.0]), SVector{2,Float64}([0.217147, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

   y = log1p(m)
   yref = MC{2}(1.3040076684760487, 1.0986122886681098, EAGO.IntervalType(0.693147, 1.60944), SVector{2,Float64}([0.30543, 0.0]), SVector{2,Float64}([0.333333, 0.0]), false)

   @test isapprox(y.cv, yref.cv; atol = mctol)
   @test isapprox(y.cc, yref.cc; atol = mctol)
   @test isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol)
   @test isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol)
   @test isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol)
   @test isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol)
   @test isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol)
   @test isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol)

end


@testset "Addition/Subtraction Constant" begin
    #EAGO.set_diff_relax(0)
    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

    X1 = X + 2.1
    X2 = 2.3 + X
    X3 = X + 2
    X4 = 3 + X
    @test +X == X
    @test X1.cv == 6.6
    @test X1.cc == 6.6
    @test X2.cv == 6.8
    @test X2.cc == 6.8
    @test X3.cv == 6.5
    @test X3.cc == 6.5
    @test X4.cv == 7.5
    @test X4.cc == 7.5

    X1 = X + Float16(2.1)
    X2 = Float16(2.3) + X
    X3 = X + Int16(2)
    X4 =  Int16(3) + X
    @test +X == X
    @test X1.cv == 6.599609375
    @test X1.cc == 6.599609375
    @test X2.cv == 6.80078125
    @test X2.cc == 6.80078125
    @test X3.cv == 6.5
    @test X3.cc == 6.5
    @test X4.cv == 7.5
    @test X4.cc == 7.5

    X1n = X - 2.1
    X2n = 2.3 - X
    X3n = X - 2
    X4n = 3 - X
    @test X1n.cv == 2.4
    @test X1n.cc == 2.4
    @test X2n.cv == -2.2
    @test X2n.cc == -2.2
    @test X3n.cv == 2.5
    @test X3n.cc == 2.5
    @test X4n.cv == -1.5
    @test X4n.cc == -1.5

    X1n = X - Float16(2.1)
    X2n = Float16(2.3) - X
    X3n = X - Int16(2)
    X4n = Int16(3) - X
    @test X1n.cv == 2.400390625
    @test X1n.cc == 2.400390625
    @test X2n.cv == -2.19921875
    @test X2n.cc == -2.19921875
    @test X3n.cv == 2.5
    @test X3n.cc == 2.5
    @test X4n.cv == -1.5
    @test X4n.cc == -1.5
end

@testset "Multiplication/Division Constant" begin

    #EAGO.set_diff_relax(0)
    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0);Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)

    X1 = X * 2.1
    X2 = 2.3 * X
    X3 = X * 2
    X4 = 3 * X
    @test X1.cv == 9.450000000000001
    @test X1.cc == 9.450000000000001
    @test X2.cv == 10.35
    @test X2.cc == 10.35
    @test X3.cv == 9.0
    @test X3.cc == 9.0
    @test X4.cv == 13.5
    @test X4.cc == 13.5

    X1 = X * Float16(2.1)
    X2 = Float16(2.3) * X
    X3 = X * Int16(2)
    X4 =  Int16(3) * X
    @test X1.cv == 9.4482421875
    @test X1.cc == 9.4482421875
    @test X2.cv == 10.353515625
    @test X2.cc == 10.353515625
    @test X3.cv == 9.0
    @test X3.cc == 9.0
    @test X4.cv == 13.5
    @test X4.cc == 13.5

    X1 = X * (-2.1)
    X2 = (-2.3) * X
    X3 = X * (-2)
    X4 = (-3) * X
    @test X1.cv == -9.450000000000001
    @test X1.cc == -9.450000000000001
    @test X2.cv == -10.35
    @test X2.cc == -10.35
    @test X3.cv == -9.0
    @test X3.cc == -9.0
    @test X4.cv == -13.5
    @test X4.cc == -13.5

    X1 = X * Float16(-2.1)
    X2 = Float16(-2.3) * X
    X3 = X * Int16(-2)
    X4 =  Int16(-3) * X
    @test X1.cv == -9.4482421875
    @test X1.cc == -9.4482421875
    @test X2.cv == -10.353515625
    @test X2.cc == -10.353515625
    @test X3.cv == -9.0
    @test X3.cc == -9.0
    @test X4.cv == -13.5
    @test X4.cc == -13.5

    a = seedg(Float64,1,2)
    b = seedg(Float64,2,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,4.0);Interval(-5.0,-3.0)])
    Y = MC{2}(-4.0,-4.0,xIBox[2],b,b,false)
    out = 1.0/Y
    out1 = 1/Y
    @test isapprox(out.cc,-0.25,atol=1E-6)
    @test isapprox(out.cv,-0.266666666,atol=1E-6)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cc_grad[2],-0.0625,atol=1E-6)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out.cv_grad[2],-0.0666667,atol=1E-6)
    @test isapprox(out.Intv.lo,-0.33333333,atol=1E-5)
    @test isapprox(out.Intv.hi,-0.199999,atol=1E-5)

    @test isapprox(out1.cc,-0.25,atol=1E-6)
    @test isapprox(out1.cv,-0.266666666,atol=1E-6)
    @test isapprox(out1.cc_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cc_grad[2],-0.0625,atol=1E-6)
    @test isapprox(out1.cv_grad[1],0.0,atol=1E-6)
    @test isapprox(out1.cv_grad[2],-0.0666667,atol=1E-6)
    @test isapprox(out1.Intv.lo,-0.33333333,atol=1E-5)
    @test isapprox(out1.Intv.hi,-0.199999,atol=1E-5)
end

@testset "Minimum/Maximum Constant" begin
    a = seed_g(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = min(3,X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,3)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = SMCg{2,Interval{Float64},Float64}(4.5,4.5,xIBox[1],a,a,false)
    out = max(5,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(X,5)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = min(3.0,X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,3.0)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(5.0,X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(X,5.0)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = min(Float16(3.0),X)
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    out = min(X,Float16(3.0))
    @test isapprox(out.cc,3.0,atol=1E-1)
    @test isapprox(out.cv,1.0909090909090908,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.545455,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,-3,atol=1E-4)
    @test isapprox(out.Intv.hi,3,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(Float16(5.0),X)
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)

    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    out = max(X,Float16(5.0))
    @test isapprox(out.cc,7.045454545454545,atol=1E-1)
    @test isapprox(out.cv,5.0,atol=1E-1)
    @test isapprox(out.cc_grad[1],0.272727,atol=1E-4)
    @test isapprox(out.cc_grad[2],0.0,atol=1E-1)
    @test isapprox(out.cv_grad[1],0.0,atol=1E-4)
    @test isapprox(out.cv_grad[2],0.0,atol=1E-1)
    @test isapprox(out.Intv.lo,5,atol=1E-4)
    @test isapprox(out.Intv.hi,8,atol=1E-4)
end

@testset "Conversion" begin
    a = seedg(Float64,1,2)
    xIBox = SVector{2,Interval{Float64}}([Interval(-3.0,8.0),Interval(-3.0,8.0)])
    X = MC{2}(4.5,4.5,xIBox[1],a,a,false)
    X1 = convert(MC{2},1)
    X2 = convert(MC{2},1.1)
    X3 = convert(MC{2},Interval(2.1,4.3))
    @test X1.cc == 1.0
    @test X1.cv == 1.0
    @test X2.cc == 1.1
    @test X2.cv == 1.1
    @test X3.cc == 4.3
    @test X3.cv == 2.1
end
