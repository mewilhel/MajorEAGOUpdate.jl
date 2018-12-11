# Relaxes nonlinear term via MidPoint Affine calculation
function MidPointAffine!(src::Optimizer,trg,n::NodeBB,r)
    ngrad = src.VariableNumber

    var = src.UpperVariables

    evaluator = src.WorkingEvaluatorBlock.evaluator
    evaluator.current_node = n
    midx = (n.UpperVar + n.LowerVar)/2.0

    println("midx: $midx")

    # Add objective
    if src.WorkingEvaluatorBlock.has_objective
        # Calculates relaxation and subgradient
        df = zeros(Float64,ngrad)
        f = MOI.eval_objective(evaluator, midx)
        MOI.eval_objective_gradient(evaluator, df, midx)

        # Add objective relaxation to model
        saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(df, var), f-sum(midx.*df))
        MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), saf)

        # Add objective cut if nonlinear
        set = MOI.LessThan(src.GlobalUpperBound)
        MOI.add_constraint(trg, saf, set)
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
            if bns.upper != Inf
                println("bns.upper: $(bns.upper)")
                constant = g[j]
                for i in 1:ngrad
                    constant -= midx[i]*dg[j,i]
                end
                saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(dg[j,:], var), 0.0)
                set = MOI.LessThan(bns.upper-constant)
                println("upper bound exists")
                println("saf: $saf, set: $set")
                MOI.add_constraint(trg, saf, set)
            end
            if bns.lower > -Inf
                println("bns.lower: $(bns.lower)")
                constant = g_cc[j]
                for i in 1:ngrad
                    constant -= midx[i]*dg_cc[j,i]
                end
                saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-dg_cc[j,:], var), 0.0)
                set = MOI.LessThan(constant)
                println("lower bound exists")
                println("saf: $saf, set: $set")
                MOI.add_constraint(trg, saf, set)
            end
        end
    end
end
