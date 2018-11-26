
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
    Reform_Imp_Y(X::Vector{Interval{Float64}},Y::Vector{Interval{Float64}},
                 pUBD::Vector{Vector{Float64}})
Reformulates the bounds on x and y, `Y`, to be bounds on `y*` for input into implicit
global optimization routine where `y* = [y_1; y_2; ... y_np,x]` returning a vector
of lower bounds, upper bounds, the state space dimension (for the opt problem), and
the entire problem dimension.
"""
function Reform_Imp_Y(xL::Vector{Float64},xU::Vector{Float64},yL::Vector{Float64},yU::Vector{Float64},P)
  nx::Int = length(xL)
  ny::Int = length(yL)
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
