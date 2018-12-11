function EAGODefault_UpperBounding!(x::Optimizer,y::NodeBB)
    if is_integer_feasible(x)
        println("upper flag 1")
        # Copies initial model into working model (if initial model isn't dummy)
        if x.InitialUpperOptimizer != DummyOptimizer()
            println("upper flag 1a")
            println("x.WorkingUpperOptimizer:")
            println("x.InitialUpperOptimizer")
            x.WorkingUpperOptimizer = deepcopy(x.InitialUpperOptimizer)
        end
        println("upper flag 2")
        # Updates variables bounds
        Update_VariableBounds_Upper!(x,y,x.WorkingUpperOptimizer)
        println("upper flag 3")
        # Optimizes the object
        #TT = stdout
        #redirect_stdout()
        MOI.optimize!(x.WorkingUpperOptimizer)
        #redirect_stdout(TT)
        println("upper flag 4")
        # Process output info and save to CurrentUpperInfo object
        termination_status = MOI.get(x.WorkingUpperOptimizer, MOI.TerminationStatus())
        println("upper flag 5")
        result_status = MOI.get(x.WorkingUpperOptimizer, MOI.PrimalStatus())
        println("upper flag 6")
        if (termination_status == MOI.Success) && (result_status == MOI.FeasiblePoint)
            println("upper flag A")
            x.CurrentUpperInfo.Feasibility = true
            x.CurrentUpperInfo.Value = MOI.get(x.WorkingUpperOptimizer, MOI.ObjectiveValue())
            x.CurrentUpperInfo.Solution[1:end] = MOI.get(x.WorkingUpperOptimizer, MOI.VariablePrimal(), x.UpperVariables)
        else
            println("upper flag B")
            x.CurrentUpperInfo.Feasibility = false
            x.CurrentUpperInfo.Value = Inf
        end
    else
        println("upper flag C")
        x.CurrentUpperInfo.Feasibility = false
        x.CurrentUpperInfo.Value = Inf
    end
end
