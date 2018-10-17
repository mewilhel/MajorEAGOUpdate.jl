function MOI.empty!(m::Optimizer)
    m = Optimizer()
end

function MOI.is_empty(m::Optimizer)
    vb = Bool[false for i=1:15]
    vb[1] = isempty(m.VariableInfo)
    vb[2] = m.StartedSolve
    vb[3] = m.OptimizationSense == MOI.FeasibilitySense
    vb[4] = m.TerminationStatusCode == MOI.OtherError
    vb[5] = m.FailedSolver == NoFailure

    vb[6] = m.LowerProblem == DummyFunction
    vb[7] = m.UpperProblem == DummyFunction
    vb[8] = m.Preprocess == DummyFunction
    vb[9] = m.Postprocess == DummyFunction
    vb[10] = m.RepeatCheck == DummyFunction
    vb[11] = m.ConvergenceCheck == DummyFunction
    vb[12] = m.TerminationCheck == DummyFunction
    vb[13] = m.NodeStorage == DummyFunction
    vb[14] = m.NodeSelection == DummyFunction
    vb[15] = m.BisectionFunction == DummyFunction
    bool_out = (sum(vb) > 0)
    return bool_out
end

function check_inbounds(m::Optimizer, vi::MOI.VariableIndex)
    num_variables = length(m.VariableInfo)
    if !(1 <= vi.value <= num_variables)
        error("Invalid variable index $vi. ($num_variables variables in the model.)")
    end
end

check_inbounds(m::Optimizer, var::MOI.SingleVariable) = check_inbounds(m, var.variable)

function check_inbounds(m::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        check_inbounds(m, term.variable_index)
    end
end

function check_inbounds(m::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        check_inbounds(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds(m, term.variable_index_1)
        check_inbounds(m, term.variable_index_2)
    end
end

function has_upper_bound(m::Optimizer, vi::MOI.VariableIndex)
    return m.VariableInfo[vi.value].has_upper_bound
end

function has_lower_bound(m::Optimizer, vi::MOI.VariableIndex)
    return m.VariableInfo[vi.value].has_lower_bound
end

function is_fixed(m::Optimizer, vi::MOI.VariableIndex)
    return m.VariableInfo[vi.value].is_fixed
end

function is_integer_feasible(m::Optimizer)
    int_feas = true
    for var in m.IntegerVar
        (0.0 < m.CurrentLowerInfo.Solution[var] < 1.0) && (int_feas = false; break)
    end
    return int_feas
end
