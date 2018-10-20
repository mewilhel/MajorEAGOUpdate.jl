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
end

# Relaxes nonlinear term via MidPoint Affine calculation
function MidPointAffine!(src::Optimizer,trg::T) where {T<:MOI.AbstractOptimizer}
end
