"""
    seed_g(T::Type,j::Int,N::Int)

Creates a `x::SVector{N,T}` object that is one at `x[j]` and zero everywhere else.
"""
seedg(T::Type,j::Int,N::Int) = SVector{N,T}([i == j ? 1.0 : 0.0 for i=1:N])

"""
    mid3(x::T,y::T,z::T)

Calculates the midpoint of three numbers returning the value and the index.
"""
function mid3(x::Float64,y::Float64,z::Float64)
  (((x>=y)&&(y>=z))||((z>=y)&&(y>=x))) && (return y,2)
  (((y>=x)&&(x>=z))||((z>=x)&&(x>=y))) && (return x,1)
  return z,3
end

"""
    mid_grad(cc_grad::SVector{N,T}, cv_grad::SVector{N,T}, id::Int64)

Takes the concave relaxation gradient 'cc_grad', the convex relaxation gradient
'cv_grad', and the index of the midpoint returned 'id' and outputs the appropriate
gradient according to McCormick relaxation rules.
"""
function mid_grad(cc_grad::SVector{N,Float64}, cv_grad::SVector{N,Float64}, id::Int) where N
  if (id == 1)
    return cc_grad
  elseif (id == 2)
    return cv_grad
  elseif (id == 3)
    return zeros(SVector{N,Float64})
  else
    error("Invalid mid3 position")
  end
end

"""
    line_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the line segment between `(x1,y1)` and `(x2,y2)` at `x = x0`.
"""
function line_seg(x0::Float64,x1::Float64,
                  y1::Float64,x2::Float64,y2::Float64)
   if (x2-x1) == zero(Float64)
     return y1
   else
     return y1*((x2-x0)/(x2-x1)) + y2*((x0-x1)/(x2-x1))
   end
end

"""
    dline_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the slope line segment between `(x1,y1)` and `(x2,y2)`
defaults to evaluating the derivative of the function if the interval is tight.
"""
function dline_seg(x0::Float64,x1::Float64,
                   y1::Float64,x2::Float64,
                   y2::Float64,d::Float64)
    if (x2 == x1)
      return d
    else
      return (y2-y1)/(x2-x1)
    end
end

"""
    grad_calc(cv::T,cc::T,int1::Int64,int2::Int64,dcv::SVector{N,T},dcc::SVector{N,T}) where {N,T}

(Sub)gradient calculation function. Takes the convex gradient, 'cv', the
concave gradient, 'cc', the mid index values 'int1,int2', and the derivative of
the convex and concave envelope functions 'dcv,dcc'.
"""
function grad_calc(cv::SVector{N,Float64},cc::SVector{N,Float64},int1::Int,int2::Int,dcv::Float64,dcc::Float64) where N
  cv_grad::SVector{N,Float64} = dcv*( int1==1 ? cv : ( int1==2 ? cv : zeros(SVector{N,Float64})))
  cc_grad::SVector{N,Float64} = dcc*( int2==1 ? cc : ( int2==2 ? cc : zeros(SVector{N,Float64})))
  return cv_grad, cc_grad
end

"""
    outer_rnd!(Intv::Interval{T})

Outer rounds the interval `Intv` by `MC_param.outer_param`.
"""
outer_rnd(Intv::IntervalType) = Intv.lo-MC_param.outer_param, Intv.hi+MC_param.outer_param

"""
    isequal(x::IntervalType,y::IntervalType,atol::Float64,rtol::Float64)

Checks that `x` and `y` are equal to with absolute tolerance `atol` and relative
tolerance `rtol`.
"""
isequal(x::IntervalType,y::IntervalType,atol::Float64,rtol::Float64) = (abs(x-y) < (atol + 0.5*abs(x+y)*rtol))

"""
    cut
"""
function cut(xL::Float64,xU::Float64,
             cv::Float64,cc::Float64,
             cv_grad::SVector{N,Float64},cc_grad::SVector{N,Float64}) where {N}
    if (cc > xU)
      cco::Float64 = xU
      cc_grado::SVector{N,Float64} = zeros(SVector{N,Float64})
    else
      cco = cc
      cc_grado = cc_grad
    end
    if (cv < xL)
      cvo::Float64 = xL
      cv_grado::SVector{N,Float64} = zeros(SVector{N,Float64})
    else
      cvo = cv
      cv_grado = cv_grad
  end
  return cvo,cco,cv_grado,cc_grado
end

function neq(x::Float64,y::Float64)
  abs(x-y)>EqualityTolerance
end
