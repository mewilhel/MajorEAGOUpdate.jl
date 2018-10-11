function MOI.add_constraint(model::Ipopt.Optimizer, v::MOI.SingleVariable, intv::MOI.Interval{Float64})
    vi = v.variable
    Ipopt.check_inbounds(model, vi)
    isnan(intv.lower) && error("Invalid lower bound $(gt.lower).")
    isnan(intv.upper) && error("Invalid lower bound $(gt.upper).")

    if Ipopt.has_lower_bound(model, vi)
        if Ipopt.has_upper_bound(model, vi)
            error("Lower & upper bound on variable $vi already exists.")
        else
            error("Lower bound on variable $vi already exists.")
        end
    end
    Ipopt.has_upper_bound(model, vi) && error("Upper bound on variable $vi already exists.")
    Ipopt.is_fixed(model, vi) && error("Variable $vi is already fixed.")

    model.variable_info[vi.value].lower_bound = intv.lower
    model.variable_info[vi.value].upper_bound = intv.upper
    model.variable_info[vi.value].has_lower_bound = true
    model.variable_info[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(vi.value)
end

function MOI.set(model::Ipopt.Optimizer, ::MOI.ConstraintSet, vi::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}, intv::MOI.Interval{Float64})
    v = vi.value
    Ipopt.check_inbounds(model, MOI.VariableIndex(v))
    isnan(intv.lower) && error("Invalid lower bound value $(intv.lower).")
    isnan(intv.upper) && error("Invalid upper bound value $(intv.upper).")

    model.variable_info[v].lower_bound = intv.lower
    model.variable_info[v].upper_bound = intv.upper
    model.variable_info[v].has_lower_bound = true
    model.variable_info[v].has_upper_bound = true
    model.variable_info[v].is_fixed = (model.variable_info[v].lower_bound == model.variable_info[v].upper_bound)

    return v
end
