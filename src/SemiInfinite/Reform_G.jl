
function Reform_HHJ(h::Function,hj::Function,xbar::Function)
  hout = (y,p) -> h(xbar,y,p)
  hjout = (y,p) -> hj(xbar,y,p)
  return hout, hjout
end


"""
--------------------------------------------------------------------------------
Function: BndProb_reform
--------------------------------------------------------------------------------
Description:
Reformates the semi-infinite constraint evaluated over the set Pset
and decision space constraints into a single  constraint function
--------------------------------------------------------------------------------
Inputs:
x:          Vector - Decision space variables
g:          Function - Decision space constraints
gSIP:       Function - Semi-infinite constraint gSIP(z,p)
Pset:       Array - Discretization set for reformulation
eps_g:      Float64 - Restriction for reformulation (for LBD, eps_g = 0)
--------------------------------------------------------------------------------
Returns:
An array corresponding to the output of the reformulated constraint.
--------------------------------------------------------------------------------
"""
function BndProb_reform(x,g::Function,gSIP::Function,Pset::Vector{Vector{Float64}},eps_g::Float64)
  temp = []
  if ~isempty(Pset)
    for i=1:length(Pset)
      push!(temp,gSIP(x,Pset[i])+eps_g)
    end
  end
  if (g != nothing)
    gval = g(x)
    for i=1:length(gval)
      push!(temp,gval[i])
    end
  end
  return temp
end


"""
    Reform_Imp_H(h::Function,x,y,pSet::Vector{Vector{Float64}},ny::Int64)
Reformulates the `h` function `h(x,y,p) = 0` into  pSet into `h(x,y*,p*) =
[h(x,y[1:ny],pSet[1]), h(x,y[(ny+1):2*ny],pSet[2]),...,h(x,y[(ny+1)*(np-1):ny*np],
pSet[end])]`for input into implicit global optimization routine.
"""
function Reform_Imp_H(h::Function,x::Vector{T},y,pUBD::Vector{Vector{Float64}},ny::Int) where T<:Real
  h_reform = zeros(T,ny*length(pUBD))
  for i=1:length(pUBD)
    h_reform[(1+ny*(i-1)):(ny*i)] = h(x,y[(1+ny*(i-1)):(ny*i)],pUBD[i])
  end
  return h_reform
end

"""
    Reform_Imp_G(h::Function,x,y,pUBD::Vector{Vector{Float64}},ny::Int64)
Reformulates the semi-infinite constraint function `g(x,y,p)` into  pSet into
`g(x,y*,p*) = [g(x,y[1:ny],pSet[1]), g(x,y[(ny+1):2*ny],pSet[2]),...,
g(x,y[(ny+1)*(np-1):ny*np],pSet[end])]` for input into implicit global
optimization routine.
"""
function Reform_Imp_G(g::Function,x::Vector{T},y,pUBD::Vector{Vector{Float64}},ny::Int,eps_g::Float64) where T<:Real
  g_reform = zeros(T,ny*length(pUBD))
  for i=1:length(pUBD)
    g_reform[(1+ny*(i-1)):(ny*i)] = g(x,y[(1+ny*(i-1)):(ny*i)],pUBD[i])+eps_g
  end
  return g_reform
end

"""
    Reform_Imp_G(h::Function,x,y,pUBD::Vector{Vector{Float64}},ny::Int64) - TO DO
Reformulates the semi-infinite constraint function `g(x,y,p)` into  pSet into
`g(x,y*,p*) = [g(x,y[1:ny],pSet[1]), g(x,y[(ny+1):2*ny],pSet[2]),...,
g(x,y[(ny+1)*(np-1):ny*np],pSet[end])]` for input into implicit global
optimization routine.
"""
function Reform_Imp_HG(h::Function,g::Function,x::Vector{T},y,pUBD::Vector{Vector{Float64}},ny::Int,gl::Int,eps_g::Float64) where T<:Real
  np = length(pUBD)
  hg_reform = zeros(T,(ny+gl)*np)
  for i=1:np
      hg_reform[(1+ny*(i-1)):(ny*i)] = h(x,y[(1+ny*(i-1)):(ny*i)],pUBD[i])+eps_g
  end
  for i=1:np
      hg_reform[(ny*np+1+gl*(i-1)):(ny*np+gl*i)] = g(x,y[(1+ny*(i-1)):(ny*i)],pUBD[i])+eps_g
  end
  return hg_reform
end

"""
    Reform_Imp_HJ(hj::Function,x,y,pUBD::Vector{Vector{Float64}},ny::Int64)
Reformulates the Jacobian w.r.t y, `hj!`, of `h(x,y,p)` into  pSet into `hj!(H,x,y*)`
for input into implicit global optimization routine where `y* = [y_1; y_2; ... y_np]`.
"""
function Reform_Imp_HJ(hj::Function,x::Vector{T},y,pUBD::Vector{Vector{Float64}},ny::Int) where T<:Real
  hj_reform = zeros(T,ny*length(pUBD),ny*length(pUBD))
  for i=1:length(pUBD)
    hj_reform[(1+ny*(i-1)):(ny*i),(1+ny*(i-1)):(ny*i)] = hj(x,y[(1+ny*(i-1)):(ny*i)],pUBD[i])
  end
  return hj_reform
end
"""
    Reform_Imp_Y(X::Vector{Interval{Float64}},Y::Vector{Interval{Float64}},
                 pUBD::Vector{Vector{Float64}})
Reformulates the bounds on x and y, `Y`, to be bounds on `y*` for input into implicit
global optimization routine where `y* = [y_1; y_2; ... y_np,x]` returning a vector
of lower bounds, upper bounds, the state space dimension (for the opt problem), and
the entire problem dimension.
"""
function Reform_Imp_Y(xL::Vector{Float64},xU::Vector{Float64},yL::Vector{Float64},yU::Vector{Float64},P)
  nx::Int = length(X)
  ny::Int = length(Y)
  np::Int = length(P)
  Y_reform_lo::Vector{Float64} = zeros(ny*np)
  Y_reform_hi::Vector{Float64} = zeros(ny*np)
  X_reform_lo::Vector{Float64} = zeros(nx)
  X_reform_hi::Vector{Float64} = zeros(nx)
  count::Int = 1
  for i=1:np
    for j=1:ny
      Y_reform_lo[count] = yL[j]
      Y_reform_hi[count] = yU[j]
      count += 1
    end
  end
  count = 1
  for j=1:nx
    X_reform_lo[count] = xL[j]
    X_reform_hi[count] = xU[j]
    count += 1
  end
  return X_reform_lo, X_reform_hi, Y_reform_lo, Y_reform_hi
end
#=
function Reform_Imp_Y(X::Vector{Interval{Float64}},Y::Vector{Interval{Float64}},P::Vector{Any})
  nx::Int64 = length(X)
  ny::Int64 = length(Y)
  Y_reform_lo::Vector{Float64} = zeros(nx)
  Y_reform_hi::Vector{Float64} = zeros(nx)
  count::Int64 = 1
  for j=1:nx
    Y_reform_lo[count] = X[j].lo
    Y_reform_hi[count] = X[j].hi
    count += 1
  end
  return Y_reform_lo, Y_reform_hi, 0, nx
end
=#
