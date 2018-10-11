
module MC_mod

using IntervalArithmetic
using StaticArrays

struct MC{N} <: Real
  cc::Float64
  cv::Float64
  cc_grad::SVector{N,Float64}
  cv_grad::SVector{N,Float64}
  Intv::IntervalType
  cnst::Bool

  function MC{N}(cc1::Float64,cv1::Float64,cc_grad1::SVector{N,Float64},cv_grad1::SVector{N,Float64},
                Intv1::IntervalType,cnst1::Bool) where {N}
    new(cc1,cv1,cc_grad1,cv_grad1,Intv1,cnst1)
  end
end

end

using StaticArrays

s = SVector{2,Float64}(1.0,2.0)
intv = MC_mod.Interval(0.0,4.0)
a = MC_mod.MC{2}(2.0,1.0,s,s,MC_mod.Interval(0.0,4.0),false)
