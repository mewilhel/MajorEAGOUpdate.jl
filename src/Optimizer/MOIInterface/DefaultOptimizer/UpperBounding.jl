function EAGODefault_UpperBounding!(x::Optimizer,y::NodeBB)
    if is_integer_feasible(x)
        # Copies initial model into working model (if initial model isn't dummy)
        if x.InitialUpperOptimizer != DummyOptimizer()
            #println("copied initial optimizer")
            #MOI.copy_to(x.WorkingUpperOptimizer,x.InitialUpperOptimizer)
            x.WorkingUpperOptimizer = deepcopy(x.InitialUpperOptimizer)
        end

        # Updates variables bounds
        Update_VariableBounds_Upper!(x,y,x.WorkingUpperOptimizer)
        println("fixed upper variables")
        x.Debug2 = x.WorkingUpperOptimizer
        # Optimizes the object
        TT = stdout
        redirect_stdout()
        MOI.optimize!(x.WorkingUpperOptimizer)
        redirect_stdout(TT)
        #redirect_stdout(tt)
        #return x.WorkingUpperOptimizer

        # Process output info and save to CurrentUpperInfo object
        termination_status = MOI.get(x.WorkingUpperOptimizer, MOI.TerminationStatus())
        #println("termination_status: $termination_status")
        objvalue = MOI.get(x.WorkingUpperOptimizer, MOI.ObjectiveValue())
        #println("objvalue: $objvalue")
        if termination_status == MOI.Success
            @assert MOI.get(x.WorkingUpperOptimizer, MOI.ResultCount()) > 0
            result_status = MOI.get(x.WorkingUpperOptimizer, MOI.PrimalStatus())
            if result_status != MOI.FeasiblePoint
                x.CurrentUpperInfo.Feasibility = false
                x.CurrentUpperInfo.Value = Inf
            end
            x.CurrentUpperInfo.Feasibility = true
            x.CurrentUpperInfo.Value = objvalue
            temp1 = MOI.get(x.WorkingUpperOptimizer, MOI.VariablePrimal(), x.UpperVariables)
            temp2 = x.CurrentUpperInfo.Solution[1:end]
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
