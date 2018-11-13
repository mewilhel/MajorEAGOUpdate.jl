function PSMCg_Kernel!(h::Function,hj::Function,z_mc::Vector{SMCg{N,V,T}},
                       aff_mc::Vector{SMCg{N,V,T}},p_mc::Vector{SMCg{N,V,T}},
                       x_mc::Vector{SMCg{N,V,T}},opt::mc_opts{T}) where {N,V<:AbstractInterval,T<:AbstractFloat}
  H,J = PreconditionSMCg(h,hj,z_mc,aff_mc,p_mc,opt)
  if (opt.CTyp == :Newton)
    if (opt.LAlg == :Dense)
      SMCg_Dense_Newton_GS!(z_mc,x_mc,J,H,opt)
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Newton_GS!(z_mc,x_mc,J,H,opt)
    end
  elseif (opt.CTyp == :Krawczyk)
    if (opt.LAlg == :Dense)
      SMCg_Dense_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
    end
  else
      error("The contractor type $(CTyp) is not currently supported. The
             contractors :Newton and :Krawczyk are currently supported.")
  end
end

function GenExpansionParams(h::Function, hj::Function,
                      X::Vector{V},
                      P::Vector{V},
                      pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat,V<:AbstractInterval}

  nxi::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = zeros(SVector{np,T})
  sone::SVector{np,T} = ones(SVector{np,T})

  x_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false) for i=1:nxi]
  xa_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false) for i=1:nxi]
  xA_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false) for i=1:nxi]
  z_mct::Vector{SMCg{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc
  z_mc = Rnd_Out_Z_All(z_mct,mc_opts.aff_correct_eps)

  p_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false) for i=1:np]
  pref_mc::Vector{SMCg{np,V,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,V,T}} = SMCg{np,V,T}[SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].Intv.lo,xA_mc[i].Intv.hi),false) for i=1:nxi]
  sto_out::Vector{Vector{SMCg{np,V,T}}} = Vector{SMCg{np,V,T}}[x_mc for j=1:(mc_opts.kmax+1)]
  sto_out[1] = copy(x_mc)
  optc = Any[szero,sone]

  for k=1:mc_opts.kmax
    PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,mc_opts)
    z_mc = Rnd_Out_Z_All(z_mc,mc_opts.aff_correct_eps)
    Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
    aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                   V(xA_mc[i].Intv.lo,xA_mc[i].Intv.hi),false) for i=1:nxi]
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  return sto_out
end

#=
function PSMCg_Kernel!(h::Function,hj::Function,z_mc::Vector{Hybrid{N,V,T}},
                       aff_mc::Vector{Hybrid{N,V,T}},p_mc::Vector{Hybrid{N,V,T}},
                       x_mc::Vector{Hybrid{N,V,T}},opt::mc_opts{T}) where {N,V<:AbstractInterval,T<:AbstractFloat}
  H,J = PreconditionSMCg(h,hj,z_mc,aff_mc,p_mc,opt)
  if (opt.CTyp == :Newton)
    if (opt.LAlg == :Dense)
      SMCg_Dense_Newton_GS!(z_mc,x_mc,J,H,opt)
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Newton_GS!(z_mc,x_mc,J,H,opt)
    end
  elseif (opt.CTyp == :Krawczyk)
    if (opt.LAlg == :Dense)
      SMCg_Dense_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
    elseif (opt.LAlg == :DenseBand)
      SMCg_DenseBand_Krawczyk_CW!(z_mc,x_mc,J,H,opt)
    end
  else
      error("The contractor type $(CTyp) is not currently supported. The
             contractors :Newton and :Krawczyk are currently supported.")
  end
end

function GenExpansionParams(h::Function, hj::Function,
                      X::Vector{V},
                      P::Vector{V},
                      pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat,V<:AbstractInterval}

  nxi::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = zeros(SVector{np,T})
  sone::SVector{np,T} = ones(SVector{np,T})

  x_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false)) for i=1:nxi]
  xa_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false)) for i=1:nxi]
  xA_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false)) for i=1:nxi]
  z_mct::Vector{HybridMC{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc
  z_mc = Rnd_Out_Z_All(z_mct,mc_opts.aff_correct_eps)

  p_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false)) for i=1:np]
  pref_mc::Vector{HybridMC{np,V,T}} = copy(p_mc)
  aff_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,V(xa_mc[i].Intv.lo,xA_mc[i].Intv.hi),false)) for i=1:nxi]
  sto_out::Vector{Vector{HybridMC{np,V,T}}} = Vector{HybridMC{np,V,T}}[x_mc for j=1:(mc_opts.kmax+1)]
  sto_out[1] = copy(x_mc)
  optc = Any[szero,sone]

  for k=1:mc_opts.kmax
    PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,mc_opts)
    z_mc = Rnd_Out_Z_All(z_mc,mc_opts.aff_correct_eps)
    Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
    aff_mc = SMCg{np,V,Float64}[SMCg{np,V,Float64}(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                   V(xA_mc[i].Intv.lo,xA_mc[i].Intv.hi),false) for i=1:nxi]
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  return sto_out
end
=#
