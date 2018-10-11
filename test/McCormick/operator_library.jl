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
