
# Defines single variable objective function
function eval_function(var::MOI.SingleVariable, x)
    return x[var.variable.value]
end

# Defines
function eval_objective(m::Optimizer, x)
    @assert !(m.nlp_data.has_objective && m.objective !== nothing)
    if m.nlp_data.has_objective
        return MOI.eval_objective(m.nlp_data.evaluator, x)
    elseif m.objective !== nothing
        return eval_function(m.objective, x)
    else
        return 0.0
    end
end
