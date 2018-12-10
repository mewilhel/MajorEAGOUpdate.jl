function Update_VariableBounds_Lower!(x::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    for i=1:x.VariableNumber
        var = x.VariableInfo[i]

        if (~var.is_integer)
            ci1,ci2,num = x.VariableIndexLow[i]
            if var.is_fixed
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.EqualTo{Float64}(y.LowerVar[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(y.UpperVar[i]))
                    MOI.set(z, MOI.ConstraintSet(), ci2, MOI.GreaterThan{Float64}(y.LowerVar[i]))
                else
                    MOI.set(z, MOI.ConstraintSet(), ci1, MOI.GreaterThan{Float64}(y.LowerVar[i]))
                end
            elseif var.has_upper_bound
                MOI.set(z, MOI.ConstraintSet(), ci1, MOI.LessThan{Float64}(y.UpperVar[i]))
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
        var_xi = MOI.SingleVariable(typevar[i])
        if var.is_integer
        else
            if var.is_fixed
                MOI.add_constraint(z, var_xi, MOI.EqualTo(y.LowerVar[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.add_constraint(z, var_xi, MOI.LessThan(y.UpperVar[i]))
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.LowerVar[i]))
                else
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.LowerVar[i]))
                end
            elseif var.has_upper_bound
                MOI.add_constraint(z, var_xi, MOI.LessThan(y.UpperVar[i]))
            end
        end
        #push!(m.VariableIndexUpp,VarTupleUpp)
    end
end

function MinusObjective!(d::JuMP.NLPEvaluator)
    if (d.has_nlobj)
        # shifts the adjacency matrix to introduce -f(x) as first element of nd array
        rowval = rowvals(d.objective.adj) .+ 1
        nzval = nonzeros(d.objective.adj)
        m, n = size(d.objective.adj)
        pushfirst!(d.objective.adj.colptr,1)
        d.objective.adj = SparseMatrixCSC{Bool,Int}(m+1,n+1,d.objective.adj.colptr,rowval,nzval)
        d.objective.adj[1,2] = true

        # shifts the node list (and parents)
        shift_nd = [JuMP.NodeData(JuMP.CALLUNIVAR,2,2)]
        for nd in d.objective.nd
            push!(shift_nd,JuMP.NodeData(nd.nodetype,nd.index,nd.parent+1))
        end
        d.objective.nd = shift_nd

        pushfirst!(d.objective.forward_storage,0.0)
        pushfirst!(d.objective.partials_storage,0.0)
        pushfirst!(d.objective.reverse_storage,0.0)
    end
end

function SetLocalNLP!(m::Optimizer)

    # Add linear and quadratic constraints to model
    for (func, set) in m.LinearLEQConstraints
         MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func, set) in m.LinearGEQConstraints
        MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func, set) in m.LinearEQConstraints
        MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func,set,ind) in m.LinearITVConstraints
        MOI.add_constraint(m.InitialUpperOptimizer, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(m.InitialUpperOptimizer, func, MOI.LessThan{Float64}(set.upper))
    end

    for (func, set) in m.QuadraticLEQConstraints
        MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func, set) in m.QuadraticGEQConstraints
        MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func, set) in m.QuadraticEQConstraints
        MOI.add_constraint(m.InitialUpperOptimizer,func,set)
    end
    for (func,set,ind) in m.QuadraticITVConstraints
        MOI.add_constraint(m.InitialUpperOptimizer, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(m.InitialUpperOptimizer, func, MOI.LessThan{Float64}(set.upper))
    end

    # Add nonlinear evaluation block
    MOI.set(m.InitialUpperOptimizer, MOI.NLPBlock(), m.NLPData)

    # Add objective sense
    #MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveSense(), m.OptimizationSense)
    MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # Specifies variables in upper problem
    m.UpperVariables = MOI.VariableIndex.(1:m.VariableNumber)

    objmult = (m.OptimizationSense == MOI.MinSense) ? 1.0 : -1.0

    # Add objective function (if any)
    if (m.Objective != nothing)
        if typeof(m.Objective) == MOI.SingleVariable
            if (m.OptimizationSense == MOI.MinSense)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.SingleVariable}(), m.Objective)
            elseif (m.OptimizationSense == MOI.MaxSense)
                neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, m.Objective.variable)], 0.0)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_var)
            else
                error("Objective sense must be MOI.MinSense or MOI.MaxSense")
            end
        elseif typeof(m.Objective) == MOI.ScalarAffineFunction{Float64}
            if (m.OptimizationSense == MOI.MinSense)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), m.Objective)
            elseif (m.OptimizationSense == MOI.MaxSense)
                neg_obj_aff_terms = []
                for term in m.Objective.terms
                    push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
                end
                neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -m.Objective.constant)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_aff)
            else
                error("Objective sense must be MOI.MinSense or MOI.MaxSense")
            end
        elseif typeof(m.Objective) == MOI.ScalarQuadraticFunction{Float64}
            if (m.OptimizationSense == MOI.MinSense)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), m.Objective)
            elseif (m.OptimizationSense == MOI.MaxSense)
                neg_obj_qda_terms = []
                for term in m.Objective.affine_terms
                    push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
                end
                neg_obj_qdq_terms = []
                for term in m.Objective.quadratic_terms
                    push!(neg_obj_aff_terms,MOI.ScalarQuadraticTerm{Float64}(-term.coefficient,term.variable_index_1,term.variable_index_2))
                end
                neg_obj_qd = ScalarQuadraticFunction{Float64}(neg_obj_qda_terms,neg_obj_qdq_terms,-m.Objective.constant)
                MOI.set(m.InitialUpperOptimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), neg_obj_qd)
            else
                error("Objective sense must be MOI.MinSense or MOI.MaxSense")
            end
        end
    else
        @assert m.NLPData != empty_nlp_data()
        if (m.OptimizationSense == MOI.MaxSense)
            MinusObjective!(m.NLPData.evaluator)
        elseif (m.OptimizationSense != MOI.MinSense)
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
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
    (m.RelaxFunction == DummyFunction)      &&   (m.RelaxFunction = RelaxModel!)
    (typeof(m.InitialRelaxedOptimizer) == DummyOptimizer) && (m.InitialRelaxedOptimizer = CPLEX.Optimizer())
    (typeof(m.WorkingRelaxedOptimizer) == DummyOptimizer) && (m.WorkingRelaxedOptimizer = CPLEX.Optimizer())
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

function Update_VariableBounds_Lower1!(m::Optimizer,y::NodeBB,z::T) where {T<:MOI.AbstractOptimizer}
    # Updates variables bounds
    typevar = m.LowerVariables
    for (i,var) in enumerate(m.VariableInfo)
        var_xi = MOI.SingleVariable(typevar[i])
        if var.is_integer
        else
            if var.is_fixed
                MOI.add_constraint(z, var_xi, MOI.EqualTo(y.LowerVar[i]))
            elseif var.has_lower_bound
                if var.has_upper_bound
                    MOI.add_constraint(z, var_xi, MOI.LessThan(y.UpperVar[i]))
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.LowerVar[i]))
                else
                    MOI.add_constraint(z, var_xi, MOI.GreaterThan(y.LowerVar[i]))
                end
            elseif var.has_upper_bound
                MOI.add_constraint(z, var_xi, MOI.LessThan(y.UpperVar[i]))
            end
        end
    end
end
