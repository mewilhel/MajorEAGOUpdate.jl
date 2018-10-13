"""
    grad(x::MC{N},j::Int) where N

sets convex and concave (sub)gradients of length `n` of `x` to be `1` at index `j`
"""
function grad(x::MC{N},j::Int) where N
  sv_grad::SVector{N,Float64} = seedg(T,j,N)
  return MC{N}(x.cc,x.cv,sv_grad,sv_grad,x.Intv,x.cnst)
end

"""
    zgrad(x::SMCg{N,T},n::Int64) where {N,T}

sets convex and concave (sub)gradients of length `n` to be zero
"""
function zgrad(x::MC{N}) where N
  grad::SVector{N,Float64} = zeros(SVector{N,Float64})
  return MC{N}(x.cc,x.cv,grad,grad,x.Intv,x.cnst)
end

"""
    Intv(x::MC)
"""
Intv(x::MC) = x.Intv
lo(x::MC) = x.Intv.lo
hi(x::MC) = x.Intv.hi
cc(x::MC) = x.cc
cv(x::MC) = x.cv
cc_grad(x::MC) = x.cc_grad
cv_grad(x::MC) = x.cv_grad
cnst(x::MC) = x.cnst


"""
    MC_Relaxations!(val::Integer)
Set differentiability of relaxations used.
"""
function MC_Differentiability!(val::Integer)
  diff_relax = val>0
  if (diff_relax>0)
    MC_param.mu = val
  else
    MC_param.mu = 0
  end
end
