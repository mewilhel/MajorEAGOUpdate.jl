# Relaxes nonlinear term via MidPoint Affine calculation
function MidPointAffine!(src::Optimizer,trg,n::NodeBB,r)
    ngrad = src.VariableNumber

    var = src.UpperVariables

    evaluator = src.WorkingEvaluatorBlock.evaluator
    evaluator.current_node = n
    midx = n.LowerVar + (n.UpperVar - n.LowerVar)/2.0

    println("src.WorkingEvaluatorBlock.has_objective: $(src.WorkingEvaluatorBlock.has_objective)")

    # Add objective
    if src.WorkingEvaluatorBlock.has_objective
        df = zeros(Float64,ngrad)
        f = MOI.eval_objective(evaluator, midx)
        println("f: $f")
        MOI.eval_objective_gradient(evaluator, df, midx)
        println("df: $df")
        MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(df, var), f-sum(midx.*df)))
    end

    # Add other affine constraints
    if length(src.WorkingEvaluatorBlock.constraint_bounds) > 0

        leng = length(src.WorkingEvaluatorBlock.constraint_bounds)
        g = zeros(Float64,leng)
        dg = zeros(Float64,leng,ngrad)
        g_cc = zeros(Float64,leng)
        dg_cc = zeros(Float64,leng,ngrad)

        MOI.eval_constraint(evaluator, g, midx)
        MOI.eval_constraint_jacobian(evaluator, dg,  midx)

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
            if bns.upper != Inf
                constant = g[j]
                for i in 1:ngrad
                    constant -= midx[i]*dg[j,i]
                end
                MOI.add_constraint(trg, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(dg[j,:], var), 0.0), MOI.LessThan(bns.upper-constant))
            end
            if bns.lower > -Inf
                constant = g_cc[j]
                for i in 1:ngrad
                    constant -= midx[i]*dg_cc[j,i]
                end
                MOI.add_constraint(trg, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-dg_cc[j,:], var), 0.0), MOI.LessThan(constant))
            end
        end
    end
end
