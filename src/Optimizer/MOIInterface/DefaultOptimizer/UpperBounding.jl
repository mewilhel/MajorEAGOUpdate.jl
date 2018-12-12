function EAGODefault_UpperBounding!(x::Optimizer,y::NodeBB)
    if is_integer_feasible(x)
        if x.UseUpperFactory
            factory = x.UpperFactory()
            x.InitialUpperOptimizer = factory
            x.UpperVariables = MOI.add_variables(x.InitialUpperOptimizer, x.VariableNumber)
            SetLocalNLP!(x)
            x.WorkingUpperOptimizer = x.InitialUpperOptimizer
        else
            if x.InitialUpperOptimizer != DummyOptimizer()
                x.WorkingUpperOptimizer = deepcopy(x.InitialUpperOptimizer)
            end
        end
        Update_VariableBounds_Upper!(x,y,x.WorkingUpperOptimizer)

        # Optimizes the object
        TT = stdout
        redirect_stdout()
        MOI.optimize!(x.WorkingUpperOptimizer)
        redirect_stdout(TT)

        # Process output info and save to CurrentUpperInfo object
        termination_status = MOI.get(x.WorkingUpperOptimizer, MOI.TerminationStatus())
        result_status = MOI.get(x.WorkingUpperOptimizer, MOI.PrimalStatus())
        if (termination_status == MOI.Success) && (result_status == MOI.FeasiblePoint)
            x.CurrentUpperInfo.Feasibility = true
            x.CurrentUpperInfo.Value = MOI.get(x.WorkingUpperOptimizer, MOI.ObjectiveValue())
            x.CurrentUpperInfo.Solution[1:end] = MOI.get(x.WorkingUpperOptimizer, MOI.VariablePrimal(), x.UpperVariables)
        else
            x.CurrentUpperInfo.Feasibility = false
            x.CurrentUpperInfo.Value = Inf
        end
    else
        x.CurrentUpperInfo.Feasibility = false
        x.CurrentUpperInfo.Value = Inf
    end
end
