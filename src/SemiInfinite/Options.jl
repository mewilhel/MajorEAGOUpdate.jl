"""SIP_opts an option type for the SIP routine. The fields it contains are as follows:
- tol: the absolute tolerance for convergence
- kmax: the maximum number of iterations
- eps_g0: initial eps_g restriction parameter, must be greater than zero
- r0: parameter used to reduce eps_g as each iteration (e.g. eps_g = eps_g/r0)
- verbosity: Sets verbosity of the solution routine. Either "None", "Normal", and "Full".
- return_hist: Sets the return type. True returns LBD/UBD and times at each iteration.
- LLP_Solve:
- LBP_Solve:
- UBP_Solve:
- LLP_Opt: LLP solver options.
- LBP_Opt: LBP solver options.
- UBP_Opt: UBP solver options.
- gSIPExp: Expression for semi-infinite inequality constraint
- hSIPExp: Expression for semi-infinite equality constraints
"""
mutable struct SIP_opts
  tol::Float64
  kmax::Int
  eps_g0::Float64
  r0::Float64
  return_hist
  hdr_intv::Int
  prnt_intv::Int
  LLP_Opt::Any
  LBP_Opt::Any
  UBP_Opt::Any
  P_LBD::Vector{Vector{Float64}}
  P_UBD::Vector{Vector{Float64}}
  Verbosity::Int
  inn_tol::Float64
end

SIP_opts() = SIP_opts(1E-3,100,1.0,1.5,false,20,1,Optimizer,Optimizer,Optimizer,
                      Vector{Float64}[],Vector{Float64}[],0,1.0E-4)
SIP_opts(X) = SIP_opts(1E-3,100,1.0,1.5,false,20,1,X,X,X,
                       Vector{Float64}[],Vector{Float64}[],0,1E-4)

#=
function set_to_default!(x::SIP_opts)
  x.LLP_Opt = EAGO_NLPSolver()
  x.LBP_Opt = EAGP_NLPSolver()
  x.UBP_Opt = EAGP_NLPSolver()
end
=#
