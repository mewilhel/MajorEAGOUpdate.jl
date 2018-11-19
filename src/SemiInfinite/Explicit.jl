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
function Explicit_SIP_Solve(f,gSIP,X_low,X_high,P_low,P_high,SIPopt::SIP_opts)

  # initializes solution
  UBDg = Inf
  LBDg = -Inf
  k = 0
  P_LBD = SIPopt.P_LBD
  P_UBD = SIPopt.P_UBD
  np = length(P_low)
  nx = length(X_low)
  pbar = P_low+(P_high - P_low)/2.0
  xbar = X_low+(X_high - X_low)/2.0
  INNg1 = Inf
  INNg2 = Inf
  feas = true
  LBP_vars = JuMP.VariableRef[]
  UBP_vars = JuMP.VariableRef[]
  LLP1_vars = JuMP.VariableRef[]
  LLP2_vars = JuMP.VariableRef[]

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
    gL_LBP = [-Inf for i=1:length(P_LBD)]
    gU_LBP = [0.0 for i=1:length(P_LBD)]
    LBP_vars, jLBP = loadscript_bnds_lbd!(SIPopt.LBP_Opt,nx, length(P_LBD), X_low, X_high,gL_LBP, gU_LBP, f, gSIP, P_LBD, LBP_vars)
    JuMP.optimize!(jLBP)

    # Process output info and save to CurrentUpperInfo object
    termination_status = JuMP.termination_status(jLBP)
    tLBP = JuMP.getsolvetime(jLBP)
    if (termination_status == MOI.Success)
        result_status = JuMP.primal_status(jLBP)
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          LBDg = JuMP.getobjectivevalue(jLBP)
          xbar = JuMP.getvalue.(LBP_vars)
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
    LLP1_vars,jLLP1 = loadscript_llp!(SIPopt.LLP_Opt, np, P_low, P_high, gSIP, xbar, MOI.MaxSense,1)
    JuMP.optimize!(jLLP1)

    # Process output info and save to CurrentUpperInfo object
    termination_status = JuMP.termination_status(jLLP1)
    tLLP = JuMP.getsolvetime(jLLP1)
    if (termination_status == MOI.Success)
        result_status = JuMP.primal_status(jLLP1)
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          INNg1 = JuMP.getobjectivevalue(jLLP1)
          pbar = JuMP.getvalue.(LLP1_vars)
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
    gL_UBP = [-Inf for i=1:length(P_UBD)]
    gU_UBP = [0.0 for i=1:length(P_UBD)]
    UBP_vars,jUBP = loadscript_bnds_ubd!(SIPopt.UBP_Opt, nx, length(P_UBD), X_low, X_high, gL_UBP, gU_UBP,  MOI.MinSense, f, gSIP, P_UBD, UBP_vars, eps_g)
    JuMP.optimize!(jUBP)

    # Process output info and save to CurrentUpperInfo object
    termination_status = JuMP.termination_status(jUBP)
    tUBP = JuMP.getsolvetime(jUBP)
    if (termination_status == MOI.Success)
        result_status = JuMP.primal_status(jUBP)
        if (result_status != MOI.FeasiblePoint)
          feas = false
        else
          feas = true
          UBD_temp = JuMP.getobjectivevalue(jUBP)
          xbar = JuMP.getvalue.(UBP_vars)
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
      LLP2_vars,jLLP2 = loadscript_llp!(SIPopt.LLP_Opt, np, P_low, P_high, gSIP, xbar, MOI.MaxSense,2)
      JuMP.optimize!(jLLP2)

      # Process output info and save to CurrentUpperInfo object
      termination_status = JuMP.termination_status(jLLP2)
      tLLP = JuMP.getsolvetime(jLLP2)
      if (termination_status == MOI.Success)
          result_status = JuMP.primal_status(jLLP2)
          if (result_status != MOI.FeasiblePoint)
            feas = false
          else
            feas = true
            INNg2 = JuMP.getobjectivevalue(jLLP2)
            pbar = JuMP.getvalue.(LLP2_vars)
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
