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

    if isa(src.Objective, MOI.SingleVariable)
        MOI.set(trg, MOI.ObjectiveFunction{MOI.SingleVariable}(), src.Objective)
    elseif isa(src.Objective, MOI.ScalarAffineFunction{Float64})
        MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), src.Objective)
    end
end
