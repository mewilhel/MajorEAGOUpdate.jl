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

# Relaxes nonlinear term via MidPoint Affine calculation
function MidPointAffine!(src::Optimizer,trg,n::NodeBB,r)
    var = src.UpperVariables

    evaluator = src.WorkingEvaluatorBlock.evaluator
    evaluator.current_node = n
    midx = n.LowerVar + (n.UpperVar - n.LowerVar)/2.0
    println("midx: $midx")

    # Add objective
    if src.WorkingEvaluatorBlock.has_objective
        df = zeros(Float64,src.VariableNumber)
        f = MOI.eval_objective(evaluator, midx)
        MOI.eval_objective_gradient(evaluator, df, midx)
        MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(df, var), f+sum(midx.*df)))
    end

    # Add other affine constraints
    if length(src.WorkingEvaluatorBlock.constraint_bounds) > 0

        leng = length(src.WorkingEvaluatorBlock.constraint_bounds)
        g = zeros(Float64,leng)
        dg = zeros(Float64,leng,src.VariableNumber)
        g_cc = zeros(Float64,leng)
        dg_cc = zeros(Float64,leng,src.VariableNumber)

        println("assigned storage")

        MOI.eval_constraint(evaluator, g, midx)
        MOI.eval_constraint_jacobian(evaluator, dg,  midx)

        println("g: $g")
        println("dg: $dg")

        # gets jacobian and gradient of convex terms
        for i in 1:length(evaluator.constraints)
            if evaluator.constraints[i].numvalued[1]
                g_cc[i] = evaluator.constraints[i].numberstorage[1]
            else
                g_cc[i] = evaluator.constraints[i].setstorage[1].cc
                dg_cc[i,:] = evaluator.constraints[i].setstorage[1].cc_grad
            end
        end

        for (j,bns) in enumerate(src.WorkingEvaluatorBlock.constraint_bounds)
            if bns.upper < Inf
                MOI.add_constraint(trg, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(dg[j,:], var), g[j]), MOI.LessThan(bns.upper))
            end
            if bns.lower > -Inf
                MOI.add_constraint(trg, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-dg_cc[j,:], var), -g_cc[j]), MOI.LessThan(-bns.lower))
            end
        end
    end
end
