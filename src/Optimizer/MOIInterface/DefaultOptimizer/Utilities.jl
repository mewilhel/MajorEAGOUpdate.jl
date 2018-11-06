function Update_VariableBounds_Lower!(x::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    for i=1:x.VariableNumber
        var = x.VariableInfo[i]
        if (~var.is_integer)
            ci1,ci2,num = x.VariableIndexLow[i]
            if var.is_fixed
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.EqualTo{Float64}(var.upper_bound))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(var.upper_bound))
                    MOI.set(z, MOI.ConstraintSet(), ci2, MOI.GreaterThan{Float64}(var.lower_bound))
                else
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.GreaterThan{Float64}(var.lower_bound))
                end
            elseif var.has_upper_bound
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(var.upper_bound))
            end
        else
            #=
            MOI.set(z, MOI.ConstraintSet(), x.VariableIndex[i], MOI.Interval(x.CurrentLowerInfo.Solution[i],
                                                                             x.CurrentLowerInfo.Solution[i]))
            if var.is_fixed
                MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.EqualTo(var.upper_bound))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.LessThan(var.upper_bound))
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.GreaterThan(var.lower_bound))
                else
                    MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.GreaterThan(var.lower_bound))
                end
            elseif var.has_upper_bound
                MOI.set(z, MOI.ConstraintSet(), MOI.VariableIndex(x.VariableNumber), MOI.LessThan(var.upper_bound)
            end
            =#
        end
    end
end

function Update_VariableBounds_Upper!(m::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    typevar = m.UpperVariables
    for (i,var) in enumerate(m.VariableInfo)
        println("i: $i")
        println("var: $var")
        var_xi = MOI.SingleVariable(typevar[i])
        println("var_xi: $var_xi")
        if var.is_integer
        else
            if var.is_fixed
                println("ran fixed")
                MOI.add_constraint(z, var_xi, MOI.EqualTo(y.LowerVar[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    println("has both bounds")
                    MOI.add_constraint(z, var_xi, MOI.LessThan(y.LowerVar[i]))
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.UpperVar[i]))
                else
                    println("has lower bound only")
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.UpperVar[i]))
                end
            elseif var.has_upper_bound
                println("has upper bound only")
                MOI.add_constraint(z, var_xi, MOI.LessThan(y.LowerVar[i]))
            end
        end
        #push!(m.VariableIndexUpp,VarTupleUpp)
    end
end

function SetLocalNLP!(m::Optimizer)
    #m.WorkingUpperOptimizer = deepcopy(m.NLPOptimizer)

    # Add variables to model
    #m.UpperVariables = MOI.add_variables(m.WorkingUpperOptimizer, m.VariableNumber)
    #for (i,var) in enumerate(m.VariableInfo)
    #    push!(m.VariableIndex, MOI.add_constraint(m.WorkingUpperOptimizer,
    #                           MOI.SingleVariable(m.UpperVariables[i]),
    #                           MOI.Interval(var.lower_bound,var.upper_bound)))
    #end

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

    # Add nonlinear evaluation block
    MOI.set(m.WorkingUpperOptimizer, MOI.NLPBlock(), m.NLPData)

    # Add objective sense
    MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveSense(), m.OptimizationSense)

    # Specifies variables in upper problem
    m.UpperVariables = MOI.VariableIndex.(1:m.VariableNumber)

    # Add objective function (if any)
    if (m.Objective != nothing)
        if typeof(m.Objective) == MOI.SingleVariable
            MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), m.Objective)
        elseif typeof(m.Objective) == MOI.ScalarAffineFunction{Float64}
            MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), m.Objective)
        elseif typeof(m.Objective) == MOI.ScalarQuadraticFunction{Float64}
            MOI.set(m.WorkingUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.Objective)
        end
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
    (typeof(m.InitialRelaxedOptimizer) == DummyOptimizer) && (m.InitialRelaxedOptimizer = Clp.Optimizer())
    (typeof(m.WorkingRelaxedOptimizer) == DummyOptimizer) && (m.WorkingRelaxedOptimizer = Clp.Optimizer())
    (typeof(m.InitialUpperOptimizer) == DummyOptimizer) && (m.InitialUpperOptimizer = Ipopt.Optimizer())
    (typeof(m.WorkingUpperOptimizer) == DummyOptimizer) && (m.WorkingUpperOptimizer = Ipopt.Optimizer())
    (typeof(m.LPOptimizer) == DummyOptimizer) && (m.LPOptimizer = Clp.Optimizer())
    (typeof(m.NLPOptimizer) == DummyOptimizer) && (m.NLPOptimizer = Ipopt.Optimizer())
end

function CreateInitialNode!(m::Optimizer)
    LowerVars = LowerBound.(m.VariableInfo)
    UpperVars = UpperBound.(m.VariableInfo)
    m.Stack[1] = NodeBB()
    m.Stack[1].LowerVar = LowerVars
    m.Stack[1].UpperVar = UpperVars
    m.MaximumNodeID += 1
end

function PushVariableBounds!(var::VariableInfo,var_xi,m)
    if var.is_integer
    else
        if var.is_fixed
            ci1 = MOI.add_constraint(m, var_xi, MOI.EqualTo(var.upper_bound))
            return ci1,ci1,1
        elseif var.has_lower_bound
            if var.has_upper_bound
                ci1 = MOI.add_constraint(m, var_xi, MOI.LessThan(var.upper_bound))
                ci2 = MOI.add_constraint(m, var_xi, MOI.GreaterThan(var.lower_bound))
                return ci1,ci2,2
            else
                ci1 = MOI.add_constraint(m, var_xi, MOI.GreaterThan(var.lower_bound))
                return ci1,ci1,1
            end
        elseif var.has_upper_bound
            ci1 = MOI.add_constraint(m, var_xi, MOI.LessThan(var.upper_bound))
            return ci1,ci1,1
        end
    end
end

function PushLowerVariables!(m::Optimizer)
    # Copies the same variables to every submodel
    MOI.add_variables(m.InitialRelaxedOptimizer, m.VariableNumber)
    x = MOI.add_variables(m.WorkingRelaxedOptimizer, m.VariableNumber)
    for (i,var) in enumerate(m.VariableInfo)
        var_xi = MOI.SingleVariable(x[i])
        PushVariableBounds!(var, var_xi, m.InitialRelaxedOptimizer)
        VarTupleLow = PushVariableBounds!(var, var_xi, m.WorkingRelaxedOptimizer)
        push!(m.VariableIndexLow,VarTupleLow)
    end
end
