# defines square operator
sqr(x::T) where T = x*x
cv_sqr_NS(x,xL,xU) = x^2
dcv_sqr_NS(x,xL,xU) = 2.0*x

cc_sqr(x,xL,xU) = (xU>xL) ? xL^2 + (xL+xU)*(x-xL) : xU^2
dcc_sqr(x,xL,xU) = (xU>xL) ? (xL+xU) : 0.0
function cv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return x^2
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && return (x^3)/xU
	return (x^3)/xL
end
function dcv_sqr(x,xL,xU)
  (0.0 <= xL || xU <= 0.0) && return 2.0*x
	((xL < 0.0) && (0.0 <= xU) && (0.0 <= x)) && (3.0*x^2)/xU
	return (3.0*x^2)/xL
end

cv_negpowpos(x,xL,xU,n) = x^n
dcv_negpowpos(x,xL,xU,n) = n*x^(n-1)
cc_negpowpos(x,xL,xU,n) = (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowpos(x,xL,xU,n) = (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

cv_pow4(x,xL,xU,n) = x^n
dcv_pow4(x,xL,xU,n) = n*x^(n-1)
cc_pow4(x,xL,xU,n) = (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_pow4(x,xL,xU,n) = (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

# convex/concave relaxation of integer powers of 1/x for negative reals
cv_negpowneg(x,xL,xU,n::Int) = isodd(n) ? cc_negpowpos(x,xL,xU,n) : x^n
dcv_negpowneg(x,xL,xU,n::Int) = isodd(n) ? dcc_negpowpos(x,xL,xU,n) : n*x^(n-1)
cc_negpowneg(x,xL,xU,n::Int) = isodd(n) ? x^n : (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowneg(x,xL,xU,n::Int) = isodd(n) ? n*x^(n-1) : (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

cv_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? cc_negpowpos(x,xL,xU,n) : x^n
dcv_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? dcc_negpowpos(x,xL,xU,n) : n*x^(n-1)
cc_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? x^n : (xU == xL) ? xU : ((xU-x)/(xU-xL))*xL^n +((x-xL)/(xU-xL))*xU^n
dcc_negpowneg(x,xL,xU,n::Float64) = isodd(Int(n)) ? n*x^(n-1) : (xU == xL) ? 0.0 : (xU^n-xL^n)/(xU-xL)

function cv_powodd(x,xL,xU,n)
  (xU <= 0.0) && return cc_pow4(x,xL,xU,n)
  (0.0 <= xL) && return x^n
  return (xL^n)*(xU-x)/(xU-xL)+(max(0.0,x))^n
end
function dcv_powodd(x,xL,xU,n)
  (xU <= 0.0) && return dcc_pow4(x,xL,xU,n)
  (0.0 <= xL) && return n*x^(n-1)
  return -(xL^n)/(xU-xL)+n*(max(0.0,x))^(n-1)
end
function cc_powodd(x,xL,xU,n)
  (xU <= 0.0) && return x^n
  (0.0 <= xL) && return cc_pow4(x,xL,xU,n)
  return (xU^n)*(x-xL)/(xU-xL)+(min(0.0,x))^n
end
function dcc_powodd(x,xL,xU,n)
  (xU <= 0.0) && return (n-1)*x^n
  (0.0 <= xL) && return dcc_pow4(x,xL,xU,n)
  return (xU^n)/(xU-xL)+n*(min(0.0,x))^(n-1)
end

function cv_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return cv_powodd(x,xL,xU,c)
      else
        return cv_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return cv_powodd(x,xL,xU,c)
        else
          return cv_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return cv_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return x^c
      (0.0 < c < 1.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c # Concave
      (c < 0.0) && return x^c # Convex
    end
  end
end

function dcv_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return dcv_powodd(x,xL,xU,c)
      else
        return dcv_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return dcv_powodd(x,xL,xU,c)
        else
          return dcv_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return dcv_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return c*x^(c-1.0)
      (0.0 < c < 1.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0) # Concave
      (c < 0.0) && return c*x^(c-1.0) # Convex
    end
  end
end

function cc_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
		println("cc pow-pos odd")
        return cc_powodd(x,xL,xU,c)
      else
		 println("cc pow-pos even")
        return cc_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
		  println("cc pow-neg dom-neg odd")
          return cc_powodd(x,xL,xU,c)
        else
		  println("cc pow-neg dom-neg even")
          return cc_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
		 println("cc pow-neg dom-pos")
        return cc_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c
      (0.0 < c < 1.0) && return x^c # Concave
      (c < 0.0) && return (xL^c)*(xU-x)/(xU-xL)+(max(0.0,x))^c # Convex
    end
  end
end

function dcc_pow(x,xL,xU,c)
  if isinteger(c)
    if (c > 0)
     if (isodd(Int(c)))
        return dcc_powodd(x,xL,xU,c)
      else
        return dcc_pow4(x,xL,xU,c)
      end
    else
      if (xU < 0.0)
        if (isodd(Int(c)))
          return dcc_powodd(x,xL,xU,c)
        else
          return dcc_negpowneg(x,xL,xU,c)
        end
      elseif (xL > 0.0)
        return dcc_negpowpos(x,xL,xU,c)
      else
        error("Function unbounded on domain")
      end
    end
  else
    if (xL >= 0.0)
      (c > 1.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0)
      (0.0 < c < 1.0) && return c*x^(c-1.0) # Concave
      (c < 0.0) && return -(xL^c)/(xU-xL)+c*(max(0.0,x))^(c-1.0) # Convex
    end
  end
end

function sqr(x::MC{N}) where N
	Intv::IntervalType = x.Intv^2
  xL::Float64 = x.Intv.lo
  xU::Float64 = x.Intv.hi
  xLc::Float64 = Intv.lo
  xUc::Float64 = Intv.hi
  eps_max::Float64 = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.hi : x.Intv.lo
	if (x.Intv.lo < zero(Float64) < x.Intv.hi)
		eps_min::Float64 = zero(Float64)
	else
		eps_min = abs(x.Intv.hi) > abs(x.Intv.lo) ?  x.Intv.lo : x.Intv.hi
	end
  midcc::Float64,cc_id::Int = mid3(x.cc,x.cv,eps_max)
  midcv::Float64,cv_id::Int = mid3(x.cc,x.cv,eps_min)
  cc::Float64 = cc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  dcc::Float64 = dcc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
	   cv::Float64 = cv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     dcv::Float64 = dcv_sqr(midcv,x.Intv.lo,x.Intv.hi)
     gdcc1::Float64 = dcc_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcv1::Float64 = dcv_sqr(x.cv,x.Intv.lo,x.Intv.hi)
     gdcc2::Float64 = dcc_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     gdcv2::Float64 = dcv_sqr(x.cc,x.Intv.lo,x.Intv.hi)
     cv_grad::SVector{N,Float64} = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
     cc_grad::SVector{N,Float64} = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
	   cv = cv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     dcv = dcv_sqr_NS(midcv,x.Intv.lo,x.Intv.hi)
     cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
     cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
     cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
  end
  return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
end

function pow_max(Intv::IntervalType,c::Float64)
  if isinteger(c)
    ci  = Int(c)
    if isodd(ci)
      return hi(Intv)
    elseif iseven(ci)
      if (abs(lo(Intv)) < abs(hi(Intv)))
        (lo(Intv) < 0.0) && (return 0.0)
        return (hi(Intv))
      else (abs(lo(Intv)) > abs(hi(Intv)))
        (hi(Intv) > 0.0) && (return 0.0)
        return (lo(Intv))
      end
    end
  end
  return hi(Intv)
end

function pow_min(Intv::IntervalType,c::Float64)
  if isinteger(c)
    ci  = Int(c)
    if isodd(ci)
      return lo(Intv)
    elseif iseven(ci)
      if (abs(lo(Intv)) < abs(hi(Intv)))
        (lo(Intv) < 0.0) && (return 0.0)
        return (lo(Intv))
      else (abs(lo(Intv)) > abs(hi(Intv)))
        (hi(Intv) > 0.0) && (return 0.0)
        return (hi(Intv))
      end
    end
  end
  return lo(Intv)
end

function pow(x::MC{N},c::Float64) where N
  if (c == 0)
	  println("ran zero arc")
    return one(x)
  elseif (c == 1)
	  println("ran one arc")
    return x
  elseif (c == 2)
	  println("ran two arc")
    return sqr(x)
  elseif (!isinteger(c) && lo(x.Intv) <= 0.0)
	  println("ran float arc")
    return exp(c*log(x))
  else
	println("ran other integer arc")
    Intv = x.Intv^c
    xL = x.Intv.lo
    xU = x.Intv.hi
    xLc = Intv.lo
    xUc = Intv.hi
    eps_max = pow_max(x.Intv,c)
    eps_min = pow_min(x.Intv,c)
    midcc,cc_id = mid3(x.cc,x.cv,eps_max)
    midcv,cv_id = mid3(x.cc,x.cv,eps_min)
    cc = cc_pow(midcc,x.Intv.lo,x.Intv.hi,c)
    cv = cv_pow(midcv,x.Intv.lo,x.Intv.hi,c)
    if (MC_param.mu >= 1)
       gdcc1 = dcc_pow(x.cv,x.Intv.lo,x.Intv.hi,c)
       gdcv1 = dcv_pow(x.cv,x.Intv.lo,x.Intv.hi,c)
       gdcc2 = dcc_pow(x.cc,x.Intv.lo,x.Intv.hi,c)
       gdcv2 = dcv_pow(x.cc,x.Intv.lo,x.Intv.hi,c)
       cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
       cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
       dcv = dcv_pow(midcv,x.Intv.lo,x.Intv.hi,c)
       dcc = dcc_pow(midcc,x.Intv.lo,x.Intv.hi,c)
       cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
       cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
       cv,cc,cv_grad,cc_grad = cut(xLc,xUc,cv,cc,cv_grad,cc_grad)
    end
    return MC{N}(cv, cc, Intv, cv_grad, cc_grad, x.cnst)
  end
end


(^)(x::MC,c::Int) = pow(x,Float64(c))
pow(x::MC,c::Int) = pow(x,Float64(c))
inv(x::MC) = pow(x,-1)
