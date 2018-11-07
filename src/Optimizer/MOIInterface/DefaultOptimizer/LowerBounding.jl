function SetDual!(x::Optimizer)
    println("started set dual")
    for (vi,VarIndxTuple) in enumerate(x.VariableIndexLow)
        (ci1,ci2,n) = VarIndxTuple
        println("VarIndxTuple: $VarIndxTuple")
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
        #MOI.copy_to(x.WorkingRelaxedOptimizer,x.InitialRelaxedOptimizer)
        x.WorkingRelaxedOptimizer = deepcopy(x.InitialRelaxedOptimizer)
        println("copied initial relaxed solver")
    end


    Update_VariableBounds_Lower!(x,y,x.WorkingRelaxedOptimizer)

    RelaxModel!(x, x.WorkingRelaxedOptimizer, y, x.Relaxation, load = false)
    println("model was relaxed")

    # Optimizes the object
    #tt = stdout
    #redirect_stdout()
    x.Debug = x.WorkingRelaxedOptimizer
    MOI.optimize!(x.WorkingRelaxedOptimizer)
    return x.WorkingRelaxedOptimizer               # CHANGE ME
    #redirect_stdout(tt)

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
        SetDual!(x)
    else
        x.CurrentLowerInfo.Feasibility = false
        x.CurrentLowerInfo.Value = Inf
    end
end
