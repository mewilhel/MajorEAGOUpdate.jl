
"""
    MC_impRelax(h::Function, hj::Function, p::Vector{SMCg{N,T}}, pmid::Vector{T},
                X::Vector{Interval{T}}, P::Vector{Interval{T}},mc_opts::mc_opts{T},
                param::Vector{Vector{SMCg{N,T}}})

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are `param` and the basic parameters of the
fixed point method are `mc_opt`.
"""
function MC_impRelax(h::Function, hj::Function, p::Vector{MC{N}}, pmid::Vector{Float64},
                     X::Vector{IntervalType}, P::Vector{IntervalType},
                     mc_opts::mc_opts,param::Vector{Vector{MC{N}}}) where N

    nx::Int = mc_opts.nx
    np::Int = mc_opts.np
    szero::SVector{np,Float64} = zeros(SVector{np,Float64})
    sone::SVector{np,Float64} = ones(SVector{np,Float64})

    x_mc::Vector{MC{np}} = [MC{np}(X[i].lo,X[i].hi,IntervalType(X[i].lo,X[i].hi),szero,szero,true) for i=1:nx]
    xa_mc::Vector{MC{np}} = [MC{np}(X[i].lo,X[i].lo,IntervalType(X[i].lo,X[i].lo),szero,szero,true) for i=1:nx]
    xA_mc::Vector{MC{np}} = [MC{np}(X[i].hi,X[i].hi,IntervalType(X[i].hi,X[i].hi),szero,szero,true) for i=1:nx]
    z_mc::Vector{MC{np}} = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc

    p_mc::Vector{MC{np}} = copy(p)
    pref_mc::Vector{MC{np}} = [MC{np}(pmid[i],pmid[i],IntervalType(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
    aff_mc::Vector{MC{np}} = [MC{np}(xa_mc[i].cv,xA_mc[i].cc,IntervalType(xa_mc[i].cv,xA_mc[i].cc),szero,szero,false) for i=1:nx]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      Affine_Exp!(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,mc_opts)
      aff_mc = MC{np}[MC{np}(x_mc[i].cv,x_mc[i].cc,IntervalType(x_mc[i].cv,x_mc[i].cc),
                             szero,szero,false) for i=1:nx]
      PMC_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    end
    return x_mc
end
