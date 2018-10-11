MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::MOI.SingleVariable)
    check_inbounds(m, func)
    m.Objective = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::MOI.ScalarAffineFunction)
    check_inbounds(m, func)
    m.Objective = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction, func::MOI.ScalarQuadraticFunction)
    check_inbounds(m, func)
    m.Objective = func
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    m.OptimizationSense = sense
    return
end

# Defines single variable objective function
function eval_function(var::MOI.SingleVariable, x)
    return x[var.variable.value]
end

# Defines affine objective function
function eval_function(aff::MOI.ScalarAffineFunction, x)
    function_value = aff.constant
    for term in aff.terms
        function_value += term.coefficient*x[term.variable_index.value]
    end
    return function_value
end

# Defines quadratic objective function
function eval_function(quad::MOI.ScalarQuadraticFunction, x)
    function_value = quad.constant
    for term in quad.affine_terms
        function_value += term.coefficient*x[term.variable_index.value]
    end
    for term in quad.quadratic_terms
        row_idx = term.variable_index_1
        col_idx = term.variable_index_2
        coefficient = term.coefficient
        if row_idx == col_idx
            function_value += 0.5*coefficient*x[row_idx.value]*x[col_idx.value]
        else
            function_value += coefficient*x[row_idx.value]*x[col_idx.value]
        end
    end
    return function_value
end

# Defines evaluation function for objective
function eval_objective(m::Optimizer, x)
    @assert !(m.NLPData.has_objective && m.Objective !== nothing)
    if m.NLPData.has_objective
        return MOI.eval_objective(m.NLPData.evaluator, x)
    elseif m.Objective !== nothing
        return eval_function(m.Objective, x)
    else
        return 0.0
    end
end
