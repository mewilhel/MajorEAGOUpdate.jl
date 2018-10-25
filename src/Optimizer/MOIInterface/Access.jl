#=
MOI.canget(::EAGOOptimizer, ::MOI.NumberOfVariables) = true
MOI.canget(::EAGOOptimizer, ::MOI.ListOfVariableIndices) = true
MOI.canget(::EAGOOptimizer, ::MOI.SolverName) = true
MOI.canget(m::EAGOOptimizer, ::MOI.TerminationStatus) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveValue) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveBound) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.RelativeGap) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.SolveTime) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.NodeCount) = m.started_solve
=#

MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [MOI.VariableIndex(i) for i in 1:length(m.VariableInfo)]
MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m.ObjectiveValue
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = length(m.VariableInfo)
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m.GlobalUpperBound
MOI.get(m::Optimizer, ::MOI.RelativeGap) = m.GlobalUpperBound
MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m.TerminationStatusCode
function MOI.get(m::Optimizer, ::MOI.SolveTime)
    IterCount = m.CurrentIterationCount
    if IterCount > 0
        return m.History.PreprocessTime[m.CurrentIterationCount] +
               m.History.PostprocessTime[m.CurrentIterationCount] +
               m.History.LowerTime[m.CurrentIterationCount] +
               m.History.UpperTime[m.CurrentIterationCount]
     else
         return 0.0
     end
 end
MOI.get(m::Optimizer, ::MOI.NodeCount) = length(m.MaximumNodeID)

##### Add primal value access (TO DO !)

# TODO: This is a bit off, because the variable primal should be available
# only after a solve. If model.inner is initialized but we haven't solved, then
# the primal values we return do not have the intended meaning.
function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.inner === nothing
        error("VariablePrimal not available.")
    end
    check_inbounds(model, vi)
    return model.inner.x[vi.value]
end

macro define_constraint_primal(function_type, set_type, prefix)
    constraint_array = Symbol(string(prefix) * "_constraints")
    offset_function = Symbol(string(prefix) * "_offset")
    quote
        function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                         ci::MOI.ConstraintIndex{$function_type, $set_type})
            if model.inner === nothing
                error("ConstraintPrimal not available.")
            end
            if !(1 <= ci.value <= length(model.$(constraint_array)))
                error("Invalid constraint index ", ci.value)
            end
            return model.inner.g[ci.value + $offset_function(model)]
        end
    end
end

@define_constraint_primal(MOI.ScalarAffineFunction{Float64},
                          MOI.LessThan{Float64}, linear_le)
@define_constraint_primal(MOI.ScalarAffineFunction{Float64},
                          MOI.GreaterThan{Float64}, linear_ge)
@define_constraint_primal(MOI.ScalarAffineFunction{Float64},
                          MOI.EqualTo{Float64}, linear_eq)
@define_constraint_primal(MOI.ScalarQuadraticFunction{Float64},
                          MOI.LessThan{Float64}, quadratic_le)
@define_constraint_primal(MOI.ScalarQuadraticFunction{Float64},
                          MOI.GreaterThan{Float64}, quadratic_ge)
@define_constraint_primal(MOI.ScalarQuadraticFunction{Float64},
                          MOI.EqualTo{Float64}, quadratic_eq)

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.LessThan{Float64}})
    if model.inner === nothing
        error("ConstraintPrimal not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !has_upper_bound(model, vi)
        error("Variable $vi has no upper bound -- ConstraintPrimal not defined.")
    end
    return model.inner.x[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.GreaterThan{Float64}})
    if model.inner === nothing
        error("ConstraintPrimal not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !has_lower_bound(model, vi)
        error("Variable $vi has no lower bound -- ConstraintPrimal not defined.")
    end
    return model.inner.x[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.EqualTo{Float64}})
    if model.inner === nothing
        error("ConstraintPrimal not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !is_fixed(model, vi)
        error("Variable $vi is not fixed -- ConstraintPrimal not defined.")
    end
    return model.inner.x[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         MOI.LessThan{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    @assert 1 <= ci.value <= length(model.linear_le_constraints)
    # TODO: Unable to find documentation in Ipopt about the signs of duals.
    # Rescaling by -1 here seems to pass the MOI tests.
    return -1 * model.inner.mult_g[ci.value + linear_le_offset(model)]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         MOI.GreaterThan{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    @assert 1 <= ci.value <= length(model.linear_ge_constraints)
    # TODO: Unable to find documentation in Ipopt about the signs of duals.
    # Rescaling by -1 here seems to pass the MOI tests.
    return -1 * model.inner.mult_g[ci.value + linear_ge_offset(model)]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         MOI.EqualTo{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    @assert 1 <= ci.value <= length(model.linear_eq_constraints)
    # TODO: Rescaling by -1 for consistency, but I don't know if this is covered
    # by tests.
    return -1 * model.inner.mult_g[ci.value + linear_eq_offset(model)]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.LessThan{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !has_upper_bound(model, vi)
        error("Variable $vi has no upper bound -- ConstraintDual not defined.")
    end
    # MOI convention is for feasible LessThan duals to be nonpositive.
    return -1 * model.inner.mult_x_U[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.GreaterThan{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !has_lower_bound(model, vi)
        error("Variable $vi has no lower bound -- ConstraintDual not defined.")
    end
    return model.inner.mult_x_L[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                         MOI.EqualTo{Float64}})
    if model.inner === nothing
        error("ConstraintDual not available.")
    end
    vi = MOI.VariableIndex(ci.value)
    check_inbounds(model, vi)
    if !is_fixed(model, vi)
        error("Variable $vi is not fixed -- ConstraintDual not defined.")
    end
    return model.inner.mult_x_L[vi.value] - model.inner.mult_x_U[vi.value]
end

function MOI.get(model::Optimizer, ::MOI.NLPBlockDual)
    if model.inner === nothing
        error("NLPBlockDual not available.")
    end
    return -1 * model.inner.mult_g[(1 + nlp_constraint_offset(model)):end]
end
