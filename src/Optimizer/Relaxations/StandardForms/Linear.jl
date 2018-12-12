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
        if (src.OptimizationSense == MOI.MinSense)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.SingleVariable}(), src.Objective)
        elseif (src.OptimizationSense == MOI.MaxSense)
            neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, src.Objective.variable)], 0.0)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_var)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    elseif isa(src.Objective, MOI.ScalarAffineFunction{Float64})
        if (src.OptimizationSense == MOI.MinSense)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), src.Objective)
        elseif (src.OptimizationSense == MOI.MaxSense)
            neg_obj_aff_terms = []
            for term in src.Objective.terms
                push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
            end
            neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -src.Objective.constant)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), neg_obj_aff)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    end
end

function ObjCutLinear!(src::Optimizer,trg) where {T<:MOI.AbstractOptimizer}
    if isa(src.Objective, MOI.SingleVariable)
        obj_set = MOI.LessThan(src.GlobalUpperBound)
        if (src.OptimizationSense == MOI.MinSense)
            MOI.add_constraint(trg, src.Objective, obj_set)
        elseif (src.OptimizationSense == MOI.MaxSense)
            neg_obj_var = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(-1.0, src.Objective.variable)], 0.0)
            MOI.add_constraint(trg, neg_obj_var, obj_set)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    elseif isa(src.Objective, MOI.ScalarAffineFunction{Float64})
        obj_set = MOI.LessThan(src.GlobalUpperBound)
        if (src.OptimizationSense == MOI.MinSense)
            MOI.add_constraint(trg, src.Objective, obj_set)
        elseif (src.OptimizationSense == MOI.MaxSense)
            neg_obj_aff_terms = []
            for term in src.Objective.terms
                push!(neg_obj_aff_terms,MOI.ScalarAffineTerm{Float64}(-term.coefficient,term.variable_index))
            end
            neg_obj_aff = MOI.ScalarAffineFunction{Float64}(neg_obj_aff_terms, -src.Objective.constant)
            MOI.add_constraint(trg, neg_obj_aff, obj_set)
        else
            error("Objective sense must be MOI.MinSense or MOI.MaxSense")
        end
    end
end
