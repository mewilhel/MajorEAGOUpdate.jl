"""
    Explicit_SIP_Solve
Solves a semi-infinite program via the algorithm presented in Mitsos2011 using
the EAGOGlobalSolver to solve the lower bounding problem, lower level problem,
and the upper bounding problem. The options for the algorithm and the global
solvers utilized are set by manipulating a SIPopt containing the options info.
Inputs:
* `f::Function`: Objective in the decision variable. Takes a single argument
                 vector that must be untyped.
* `gSIP::Function`: The semi-infinite constraint. Takes two arguments: the first
                    being a vector containing the decision variable and the
                    second being a vector containing the uncertainity
                    variables. The function must be untyped.
* `X::Vector{Interval}`: Box constraints for decision variables
* `P::Vector{Interval}`: Box constraints for uncertainty variables
* `SIPopt::SIP_opts`: Option type containing problem information
Returns: A SIP_result composite type containing solution information.
"""
function Explicit_SIP_Solve(f,gSIP,X,P,SIPopt::SIP_opts)

  # initializes solution
  UBDg = Inf
  LBDg = -Inf
  k = 0
  P_LBD = SIPopt.P_LBD
  P_UBD = SIPopt.P_UBD
  np = length(P)
  nx = length(X)
  P_low = [P[i].lo for i=1:np]
  P_high = [P[i].hi for i=1:np]
  X_low = [X[i].lo for i=1:nx]
  X_high = [X[i].hi for i=1:nx]
  pbar = mid.(P)
  xbar = mid.(X)
  INNg1 = Inf
  INNg2 = Inf
  feas = true

  sip_sto = SIP_result()

  # checks inputs
  if (SIPopt.r0<=1)
    error("r0 must be greater than 1")
  elseif (SIPopt.eps_g0<=0)
    error("eps_g must be greater than 0")
  else
    eps_g = SIPopt.eps_g0
    r = SIPopt.r0
  end

  ##### checks for convergence #####
  for k=1:SIPopt.kmax

    ##### check for termination #####
    if (abs(UBDg-LBDg)<SIPopt.tol)
      println("Algorithm Converged")
      break
    end

    ##### lower bounding problem #####
    gLBP = x -> BndProb_reform(x,nothing,gSIP,P_LBD,0.0) # reformulate constraints
    gL_LBP = [-Inf for i=1:length(P_LBD)]
    gU_LBP = [0.0 for i=1:length(P_LBD)]
    mLBP = deepcopy(SIPopt.LBP_Opt)
    LBP_vars = EAGO.loadscript!(mLBP, nx, length(P_LBD), X_low, X_high,gL_LBP, gU_LBP, MOI.MinSense, f, gLBP)
    MOI.optimize!(mLBP)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(mLBP, MOI.TerminationStatus())
    tLBP = MOI.get(mLBP, MOI.SolveTime())
    if (termination_status == MOI.Success)
        @assert MOI.get(mLBP, MOI.ResultCount()) > 0
        result_status = MOI.get(mLBP, MOI.PrimalStatus())
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          LBDg = MOI.get(mLBP, MOI.ObjectiveValue())
          xbar = MOI.get(mLBP, MOI.VariablePrimal(), LBP_vars)
        end
    else
      error("Optimizer did not successfully terminate")
    end

    sip_sto.LBP_time += tLBP
    sip_sto.LBD = LBDg
    sip_sto.xbar = xbar

    if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
      println("solved LBD: ",LBDg," ",xbar," ",feas)
    end
    if (~feas)
      println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
      sip_sto.feas = false
      return sip_sto
    end

    ##### inner program #####
    mLLP1 = deepcopy(SIPopt.LLP_Opt)
    LLP1_vars = EAGO.loadscript!(mLLP1, np, 0, P_low, P_high,
                              Float64[], Float64[], MOI.MinSense, p -> -gSIP(xbar,p), [])
    MOI.optimize!(mLLP1)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(mLLP1, MOI.TerminationStatus())
    tLLP = MOI.get(mLLP1, MOI.SolveTime())
    if (termination_status == MOI.Success)
        @assert MOI.get(mLLP1, MOI.ResultCount()) > 0
        result_status = MOI.get(mLLP1, MOI.PrimalStatus())
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          INNg1 = MOI.get(mLLP1, MOI.ObjectiveValue())
          pbar = MOI.get(mLLP1, MOI.VariablePrimal(), LLP1_vars)
        end
    else
      error("Optimizer did not successfully terminate")
    end

    sip_sto.LBP_time += tLBP
    sip_sto.LBD = LBDg
    sip_sto.xbar = xbar

    INNg1 = - INNg1
    sip_sto.LLP_time += tLLP
    if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
      println("solved INN #1: ",INNg1," ",pbar," ",feas)
    end
    if (INNg1+SIPopt.inn_tol<=0)
      xstar = xbar
      UBDg = LBDg
      return LBDg,UBDg,sip_sto
    else
      push!(P_LBD,pbar)
    end

    ##### upper bounding problem #####
    gUBP = x -> BndProb_reform(x,nothing,gSIP,P_UBD,eps_g)
    gL_UBP = [-Inf for i=1:length(P_UBD)]
    gU_UBP = [0.0 for i=1:length(P_UBD)]

    mUBP = deepcopy(SIPopt.UBP_Opt)
    mUBP_vars = EAGO.loadscript!(mUBP, nx, length(P_UBD), X_low, X_high,
                                 gL_UBP, gU_UBP, MOI.MinSense, f, gUBP)
    MOI.optimize!(mUBP)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(mUBP, MOI.TerminationStatus())
    tUBP = MOI.get(mUBP, MOI.SolveTime())
    if (termination_status == MOI.Success)
        @assert MOI.get(mUBP, MOI.ResultCount()) > 0
        result_status = MOI.get(mUBP, MOI.PrimalStatus())
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          UBD_temp = MOI.get(mUBP, MOI.ObjectiveValue())
          xbar = MOI.get(mUBP, MOI.VariablePrimal(), mUBP_vars)
        end
    else
      error("Optimizer did not successfully terminate")
    end

    sip_sto.UBP_time += tUBP
    sip_sto.UBD = UBD_temp
    sip_sto.xbar = xbar

    if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
      println("solved UBD: ",UBD_temp," ",xbar," ",feas)
    end
    if (feas)

      ##### inner program #####
      mLLP2 = deepcopy(SIPopt.LLP_Opt)
      mLLP2_vars = EAGO.loadscript!(mLLP2, np, 0, P_low, P_high,
                                   Float64[], Float64[], MOI.MinSense, p -> -gSIP(xbar,p), [])
      MOI.optimize!(mLLP2)

      # Process output info and save to CurrentUpperInfo object
      termination_status = MOI.get(mLLP2, MOI.TerminationStatus())
      tLLP = MOI.get(mLLP2, MOI.SolveTime())
      if (termination_status == MOI.Success)
          @assert MOI.get(mLLP2, MOI.ResultCount()) > 0
          result_status = MOI.get(mLLP2, MOI.PrimalStatus())
          if (result_status != MOI.FeasiblePoint)
            feas = false
          else
            feas = true
            INNg2 = MOI.get(mLLP2, MOI.ObjectiveValue())
            pbar = MOI.get(mLLP2, MOI.VariablePrimal(), mLLP2_vars)
          end
      else
        error("Optimizer did not successfully terminate")
      end

      sip_sto.LLP_time += tLLP
      INNg2 = - INNg2
      if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
        println("solved INN #2: ",INNg2," ",pbar," ",feas)
      end
      if (INNg2+SIPopt.inn_tol<0)
        if (UBD_temp <= UBDg)
          UBDg = UBD_temp
          xstar = xbar
        end
        eps_g = eps_g/r
      else
        push!(P_UBD,pbar)
      end
    else
      eps_g = eps_g/r
    end

    print_int!(SIPopt,k,LBDg,UBDg,eps_g,r)
    sip_sto.k += k
  end

  return sip_sto
end
