function PreconditionMC(h::Function,
                          hj::Function,
                          z_mc::Vector{MC{N}},
                          aff_mc::Vector{MC{N}},
                          p_mc::Vector{MC{N}},
                          opt::mc_opts) where {N}
    if (opt.LAlg == :DenseBand)
        H,J = MC_DenseBand_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    elseif (opt.LAlg == :DenseBlockDiag)
        H,J = MC_DenseBlockDiag_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    elseif (opt.LAlg == :Dense)
        H,J = MC_Dense_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    else
        error("Unsupported Linear Algebra Style")
    end
    return H,J
end

function MC_Dense_Precondition!(h::Function,
                             hj::Function,
                             z_mc::Vector{MC{N}},
                             aff_mc::Vector{MC{N}},
                             p_mc::Vector{MC{N}},
                             opt::mc_opts) where N
    #println("opt.nx: $(opt.nx)")
    H::Vector{MC{N}} = h(z_mc,p_mc)
    #println("size H: $(size(H))")
    J::VecOrMat{MC{N}} = hj(aff_mc,p_mc)
    #println("size J: $(size(J))")
    Y = [mid(J[i,j].Intv) for i=1:opt.nx, j=1:opt.nx]
    #println("size Y: $(size(Y))")
    if (opt.nx == 1)
        YH::Vector{MC{N}} = H/Y[1]
        YJ::VecOrMat{MC{N}} = J/Y[1]
    else
        F = lufact(Y)
        YH = F\H
        YJ = F\J
    end
    return YH,YJ
end


"""
    Smooth_Cut(x_mc::MC{N},x_mc_int::MC{N})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable fashion.
"""
function Smooth_Cut(x_mc::MC{N},x_mc_int::MC{N}) where N
  t_cv::MC{N} = x_mc + max(0.0,x_mc_int-x_mc)
  t_cc::MC{N} = x_mc + min(0.0,x_mc-x_mc_int)
  return MC{N}(t_cv.cv,t_cc.cc,(x_mc.Intv ∩ x_mc_int.Intv),
               t_cv.cv_grad,t_cc.cc_grad,(t_cv.cnst && t_cc.cnst))
end

"""
    Final_Cut(x_mc::MC{N},x_mc_int::MC{N})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable or nonsmooth fashion as specified by the `MC_param.mu flag`.
"""
function Final_Cut(x_mc::MC{N},x_mc_int::MC{N}) where N
  if (MC_param.mu < 1)
    Intv = x_mc.Intv ∩ x_mc_int.Intv
    if (x_mc.cc < x_mc_int.cc)
      cc = x_mc.cc
      cc_grad::SVector{N,Float64} = x_mc.cc_grad
    else
      #println("assign cc_int: $(x_mc_int.cc)")
      cc = x_mc_int.cc
      cc_grad = x_mc_int.cc_grad
    end
    #println("cc: $cc")
    if (x_mc.cv > x_mc_int.cv)
      cv = x_mc.cv
      cv_grad::SVector{N,Float64} = x_mc.cv_grad
    else
      cv = x_mc_int.cv
      cv_grad = x_mc_int.cv_grad
    end
    #println("cv: $cv")
    x_mc::MC{N} = MC{N}(cv,cc,(x_mc.Intv ∩ x_mc_int.Intv),cv_grad,cc_grad,x_mc.cnst)
  else
    x_mc = Smooth_Cut(x_mc,x_mc_int)
  end
  return x_mc
end

"""
    Rnd_Out_Z_Intv(z_mct::MC{N},epsvi::Float64)

Rounds the interval of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_Intv(z_mct::Vector{MC{N}},epsvi::Float64) where N
  return [MC{N}(z_mct[i].cv, z_mct[i].cc, IntervalType(z_mct[i].Intv.lo-epsvi, z_mct[i].Intv.hi+epsvi),
                z_mct[i].cv_grad, z_mct[i].cc_grad, z_mct[i].cnst) for i=1:length(z_mct)]
end

"""
    Rnd_Out_Z_All(z_mct::Vector{MC{N}},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_All(z_mct::Vector{MC{N}},epsvi::Float64) where {N}
  return [MC{N}(z_mct[i].cv-epsvi, z_mct[i].cc+epsvi, IntervalType(z_mct[i].Intv.lo-epsvi, z_mct[i].Intv.hi+epsvi),
                z_mct[i].cv_grad, z_mct[i].cc_grad, z_mct[i].cnst) for i=1:length(z_mct)]
end

"""
    Rnd_Out_H_Intv(z_mct::Vector{MC{N}},Y_mct::Array{MC{N},2},epsvi::S)

Rounds the interval bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_All(z_mct::Vector{MC{N}},Y_mct::Array{MC{N},2},epsvi::Float64) where N
  temp1::Vector{MC{N}} = [MC{N}(z_mct[i].cv-epsvi, z_mct[i].cc+epsvi, IntervalType(z_mct[i].Intv.lo-epsvi, z_mct[i].Intv.hi+epsvi),
                                z_mct[i].cv_grad, z_mct[i].cc_grad, z_mct[i].cnst) for i=1:length(z_mct)]
  temp2::Array{MC{N},2} = [MC{N}(Y_mct[i,j].cv-epsvi, Y_mct[i,j].cc+epsvi,
                                 Y_mct[i,j].cv_grad, Y_mct[i,j].cc_grad, IntervalType(Y_mct[i,j].Intv.lo-epsvi, Y_mct[i,j].Intv.hi+epsvi),
                                 Y_mct[i,j].cnst) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

"""
    Rnd_Out_H_All(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_Intv(z_mct::Vector{MC{N}},Y_mct::Array{MC{N},2},epsvi::Float64) where N
  temp1::Vector{MC{N}} = [MC{N}(z_mct[i].cv, z_mct[i].cc, IntervalType(z_mct[i].Intv.lo-epsvi, z_mct[i].Intv.hi+epsvi),
                                z_mct[i].cv_grad, z_mct[i].cc_grad, z_mct[i].cnst) for i=1:length(z_mct)]
  temp2::Array{MC{N},2} = [MC{N}(Y_mct[i,j].cv, Y_mct[i,j].cc, IntervalType(Y_mct[i,j].Intv.lo-epsvi, Y_mct[i,j].Intv.hi+epsvi),
                                 Y_mct[i,j].cv_grad, Y_mct[i,j].cc_grad, Y_mct[i,j].cnst) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end



#=
"""
    Precondition(hm::Vector{SMCg{N,T}},hJm::Union{Vector{SMCg{N,T}},Array{SMCg{N,T},2}},
                 Y::Union{Vector{T},Array{T,2}},nx::Int64)

Preconditions `hm` and `hJm` by `Y` in place where all dimensions are `nx`.
"""
function Precondition!(hm::Vector{MC{N}},hJm::Vector{MC{N}},
                      Y::Vector{T},nx::Int64) where N
  S1::MC{N},S2::MC{N} = zero(MC{N}),zero(MC{N})
  for i=1:nx
    S2 = zero(MC{N})
    for j=1:nx
      S1 = zero(MC{N})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{MC{N}},hJm::Vector{MC{N}},
                      Y::Array{T,2},nx::Int64) where N
  S1::MC{N},S2::MC{N} = zero(MC{N}),zero(MC{N})
  for i=1:nx
    S2 = zero(MC{N})
    for j=1:nx
      S1 = zero(MC{N})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{MC{N}},hJm::Array{MC{N},2},
                      Y::Array{T,2},nx::Int64) where N
  S1::MC{N},S2::MC{N} = zero(MC{N}),zero(MC{N})
  for i=1:nx
    S2 = zero(MC{N})
    for j=1:nx
      S1 = zero(MC{N})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end

function Precondition!(hm::Vector{MC{N}},hJm::Array{MC{N},2},
                      Y::Vector{T},nx::Int64) where N
  S1::MC{N},S2::MC{N} = zero(MC{N}),zero(MC{N})
  for i=1:nx
    S2 = zero(MC{N})
    for j=1:nx
      S1 = zero(MC{N})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end
=#
