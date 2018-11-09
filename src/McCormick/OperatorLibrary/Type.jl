#abstract type AbstractMC <: Real end
"""
    MC{N}

`MC` is the smooth McCormick (w/ gradient) structure which is used to overload
standard calculations. The fields are:
* `cc::Float64`: Concave relaxation
* `cv::Float64`: Convex relaxation
* `cc_grad::SVector{N,Float64}`: (Sub)gradient of concave relaxation
* `cv_grad::SVector{N,Float64}`: (Sub)gradient of convex relaxation
* `Intv::V`: Interval bounds
* `cnst::Bool`: Flag for whether the bounds are constant
"""
struct MC{N} <: Real
  cv::Float64
  cc::Float64
  Intv::IntervalType
  cv_grad::SVector{N,Float64}
  cc_grad::SVector{N,Float64}
  cnst::Bool

  function MC{N}(cv1::Float64,cc1::Float64,Intv1::IntervalType,
                 cv_grad1::SVector{N,Float64},cc_grad1::SVector{N,Float64},cnst1::Bool) where {N}
    new(cv1,cc1,Intv1,cv_grad1,cc_grad1,cnst1)
  end
end

"""MC(y::Interval) initializes the differentiable McCormick object with an interval
"""
MC{N}(y::IntervalType) where N = MC{N}(y.hi,y.lo,y,SVector{N,Float64}(zeros(Float64,N)),SVector{N,Float64}(zeros(Float64,N)),true)
MC{N}(val,Intv::IntervalType) where N = MC{N}(val,val,Intv,SVector{N,Float64}(zeros(Float64,N)),SVector{N,Float64}(zeros(Float64,N)),true)
