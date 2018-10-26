function EAGODefault_LowerBounding!(x::Optimizer,y::NodeBB)

    # Copies initial model into working model (if initial model isn't dummy)
    # A dummy model is left iff all terms are relaxed
    if x.InitialRelaxedOptimizer != DummyOptimizer()
        #MOI.copy_to(x.WorkingRelaxedOptimizer,x.InitialRelaxedOptimizer)
        x.WorkingRelaxedOptimizer = deepcopy(x.InitialRelaxedOptimizer)
    end

    Update_VariableBounds_Lower!(x,y,x.WorkingRelaxedOptimizer)

    RelaxModel!(x, x.WorkingRelaxedOptimizer, y, x.Relaxation, load = false)

    MOI.optimize!(x.WorkingRelaxedOptimizer)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(x.WorkingRelaxedOptimizer, MOI.TerminationStatus())
    objvalue = MOI.get(x.WorkingRelaxedOptimizer, MOI.ObjectiveValue())
    if (termination_status == MOI.Success)
        @assert MOI.get(x.WorkingRelaxedOptimizer, MOI.ResultCount()) > 0
        result_status = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
        if (result_status != MOI.FeasiblePoint)
            x.CurrentLowerInfo.Feasibility = false
            x.CurrentLowerInfo.Value = Inf
        end
        x.CurrentLowerInfo.Feasibility = true
        x.CurrentLowerInfo.Value = objvalue
        x.CurrentLowerInfo.Solution[1:end] = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(), x.LowerVariables)
        for (vi,ci) in x.VariableIndex
            tc = typeof(ci)
            if (tc == MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}})
                x.CurrentLowerInfo.LowerVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci)
            elseif (tc == MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}})
                x.CurrentLowerInfo.UpperVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci)
            elseif (tc == MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{Float64}})
                x.CurrentLowerInfo.LowerVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci)
                x.CurrentLowerInfo.UpperVarDual[vi] = MOI.get(x.WorkingRelaxedOptimizer, MOI.ConstraintDual(), ci)
            else
                error("Not a variable constraint index.")
            end
        end
    else
        x.CurrentLowerInfo.Feasibility = false
        x.CurrentLowerInfo.Value = Inf
    end
end
