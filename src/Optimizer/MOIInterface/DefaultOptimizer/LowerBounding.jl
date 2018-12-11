function SetDual!(x::Optimizer)
    for (vi,VarIndxTuple) in enumerate(x.VariableIndexLow)
        (ci1,ci2,n) = VarIndxTuple
        if (n == 2)
            if isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}})
                x.CurrentLowerInfo.LowerVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci1)
                x.CurrentLowerInfo.UpperVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci2)
            else
                x.CurrentLowerInfo.LowerVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci2)
                x.CurrentLowerInfo.UpperVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci1)
            end
        else
            if isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}})
                x.CurrentLowerInfo.LowerVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci1)
            elseif isa(ci1,MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}})
                x.CurrentLowerInfo.UpperVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci1)
            else
                error("Not a variable constraint index.")
            end
        end
    end
end

function EAGODefault_LowerBounding!(x::Optimizer,y::NodeBB)

    # Copies initial model into working model (if initial model isn't dummy)
    # A dummy model is left iff all terms are relaxed
    if x.InitialRelaxedOptimizer != DummyOptimizer()
        x.WorkingRelaxedOptimizer = deepcopy(x.InitialRelaxedOptimizer)
    end

    Update_VariableBounds_Lower1!(x,y,x.WorkingRelaxedOptimizer)

    x.RelaxFunction(x, x.WorkingRelaxedOptimizer, y, x.Relaxation, load = true)
    x.RelaxFunction(x, x.WorkingRelaxedOptimizer, y, x.Relaxation, load = false)

    # Optimizes the object
    #tt = stdout
    #redirect_stdout()
    MOI.optimize!(x.WorkingRelaxedOptimizer)
    #redirect_stdout(tt)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(x.WorkingRelaxedOptimizer, MOI.TerminationStatus())
    if (termination_status == MOI.Success)
        result_status = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
        if (result_status == MOI.FeasiblePoint)
            x.CurrentLowerInfo.Feasibility = true
            x.CurrentLowerInfo.Value = MOI.get(x.WorkingRelaxedOptimizer, MOI.ObjectiveValue())
            x.CurrentLowerInfo.Solution[1:end] = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(), x.LowerVariables)
            SetDual!(x)
        elseif (result_status == MOI.InfeasiblePoint) || (result_status == MOI.InfeasibilityCertificate)
            x.CurrentLowerInfo.Feasibility = false
            x.CurrentLowerInfo.Value = -Inf
        end
    elseif (termination_status == MOI.InfeasibleNoResult)
        x.CurrentLowerInfo.Feasibility = false
        x.CurrentLowerInfo.Value = -Inf
    else
        error("Lower problem returned ResultStatusCode = $(result_status). Lower problem must return a
              MOI.FeasiblePoint, MOI.InfeasiblePoint, or MOI.InfeasibilityCertificate status code.")
    end
end
