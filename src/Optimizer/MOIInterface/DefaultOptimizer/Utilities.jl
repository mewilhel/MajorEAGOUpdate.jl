function Update_VariableBounds_Lower!(x::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    for i=1:x.VariableNumber
        MOI.set(z, MOI.ConstraintSet(), x.VariableIndex[i], MOI.Interval(y.LowerVar[i],y.UpperVar[i]))
    end
end

function Update_VariableBounds_Upper!(x::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    for i=1:x.VariableNumber
        if (~x.VariableInfo[i].is_integer)
            MOI.set(z, MOI.ConstraintSet(), x.VariableIndex[i], MOI.Interval(y.LowerVar[i],y.UpperVar[i]))
        else
            MOI.set(z, MOI.ConstraintSet(), x.VariableIndex[i], MOI.Interval(x.CurrentLowerInfo.Solution[i],
                                                                             x.CurrentLowerInfo.Solution[i]))
        end
    end
end

function SetLocalNLP!(m::Optimizer)
    m.WorkingUpperOptimizer = deepcopy(m.NLPOptimizer)

    # Add variables to model
    m.UpperVariables = MOI.add_variables(m.WorkingUpperOptimizer, m.VariableNumber)
    for (i,var) in enumerate(m.VariableInfo)
        push!(m.VariableIndex, MOI.add_constraint(m.WorkingUpperOptimizer,
                               MOI.SingleVariable(m.UpperVariables[i]),
                               MOI.Interval(var.lower_bound,var.upper_bound)))
    end

    # Add linear and quadratic constraints to model
    for (func, set) in m.LinearLEQConstraints
         MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for (func, set) in m.LinearGEQConstraints
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for (func, set) in m.LinearEQConstraints
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for (func, set) in m.QuadraticLEQConstraints
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for (func, set) in m.QuadraticGEQConstraints
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for (func, set) in m.QuadraticEQConstraints
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end
    for bound in m.NLPData.constraint_bounds
        MOI.add_constraint(m.WorkingUpperOptimizer,func,set)
    end

    # Add nonlinear evaluation block
    MOI.set(m.WorkingUpperOptimizer, MOI.NLPBlock(), m.NLPData)

    # Add objective sense
    MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveSense(), m.OptimizationSense)

    # Add objective function (if any)
    if typeof(m.Objective) == MOI.SingleVariable
        MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), m.Objective)
    elseif typeof(m.Objective) == MOI.ScalarAffineFunction{Float64}
        MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), m.Objective)
    elseif typeof(m.Objective) == MOI.ScalarQuadraticFunction{Float64}
        MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.Objective)
    else
        @assert m.NLPData != empty_nlp_data()
    end
end

# Sets functions to default if not otherwise specified
function SetToDefault!(m::Optimizer)
    (m.LowerProblem == DummyFunction)       &&   (m.LowerProblem = EAGODefault_LowerBounding!)
    (m.UpperProblem == DummyFunction)       &&   (m.UpperProblem = EAGODefault_UpperBounding!)
    (m.Preprocess == DummyFunction)         &&   (m.Preprocess = EAGODefault_PreProcess!)
    (m.Postprocess == DummyFunction)        &&   (m.Postprocess = EAGODefault_PostProcess!)
    (m.RepeatCheck == DummyFunction)        &&   (m.RepeatCheck = DefaultRepeatCheck)
    (m.ConvergenceCheck == DummyFunction)   &&   (m.ConvergenceCheck = DefaultConvergenceCheck)
    (m.TerminationCheck == DummyFunction)   &&   (m.TerminationCheck = DefaultTerminationCheck)
    (m.NodeStorage == DummyFunction)        &&   (m.NodeStorage = DefaultStorage)
    (m.NodeSelection == DummyFunction)      &&   (m.NodeSelection = NodeSelectBest!)
    (m.BisectionFunction == DummyFunction)  &&   (m.BisectionFunction = ContinuousRelativeBisect)
    (m.CutCondition == DummyFunction)       &&   (m.CutCondition = DefaultCutCondition)
    (m.AddCut == DummyFunction)             &&   (m.AddCut = DefaultAddCut)
end

function CreateInitialNode!(m::Optimizer)
    LowerVars = LowerBound.(m.VariableInfo)
    UpperVars = UpperBound.(m.VariableInfo)
    m.Stack[1] = NodeData(LowerVars,UpperVars,-Inf,Inf,0)
    m.MaximumNodeID += 1
end
