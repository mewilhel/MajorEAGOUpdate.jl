function RelaxLinear!(src::Optimizer,trg::T) where {T<:MOI.AbstractOptimizer}
    for (func,set,ind) in src.LinearLEQConstraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.LinearGEQConstraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.LinearEQConstraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.LinearITVConstraints
        MOI.add_constraint(trg, func, MOI.GreaterThan{Float64}(set.lower))
        MOI.add_constraint(trg, func, MOI.LessThan{Float64}(set.upper))
    end

    if isa(src.Objective, MOI.SingleVariable)
        MOI.set(trg, MOI.ObjectiveFunction{MOI.SingleVariable}(), src.Objective)
    elseif isa(src.Objective, MOI.ScalarAffineFunction{Float64})
        MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), src.Objective)
    end
end
