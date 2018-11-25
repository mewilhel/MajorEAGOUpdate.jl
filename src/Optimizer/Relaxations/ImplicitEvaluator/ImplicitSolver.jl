# Modifies functions post initial relaxation to use appropriate nlp evaluator
function ImplicitMod(opt::Optimizer,args)
    ImpLowerEval = args[1]; ImpUpperEval = args[2]
    src.WorkingEvaluatorBlock = MOI.NLPBlockData(opt.NLPData.constraint_bounds,
                                                 ImpLowerEval,
                                                 opt.NLPData.has_objective)
    # copy working evaluator into block if nonlinear block is needed
    if MOI.supports(opt, MOI.NLPBlock())
        if ~isempty(opt.NonlinearVariable)
            opt.InitialRelaxedOptimizer.nlp_data = opt.WorkingEvaluatorBlock
        end
    end
end

function SolveImplicit(h::Function,
                       xl::Vector{Float64}, xu::Vector{Float64},
                       pl::Vector{Float64}, pu::Vector{Float64},
                       opt::Optimizer; f::Function = f, g::Function = g, hj::Function = hj,
                       user_sparsity::Vector{Tuple{Int64,Int64}} = sparse_pattern)

    #
    @assert length(pl) == length(pu)
    @assert length(xl) == length(xu)
    np = length(pl); nx = length(xl); ng = length(g(pl,xl))

    # Sets most routines to default
    SetToDefault!(opt)
    opt.BisectionFunction = ImplicitBisection

    # add variables to lower, upper, and EAGO models
    var_EAGO = MOI.add_variables(opt, np+nx)
    for i in 1:np
        MOI.add_constraint(opt, var_EAGO[i], MOI.LessThan(pl[i]))
        MOI.add_constraint(opt, var_EAGO[i], MOI.GreaterThan(pu[i]))
    end
    for j in 1:nx
        MOI.add_constraint(opt, var_EAGO[j+np], MOI.LessThan(xl[j]))
        MOI.add_constraint(opt, var_EAGO[j+np], MOI.GreaterThan(xu[j]))
    end

    # Build the lower evaluator
    ImpLowerEval = ImplicitLowerEvaluator()
    build_lower_evaluator!(ImpLowerEval, obj = f, constr = g, impfun = h,
                           nx = nx, np = np, user_sparse = sparse_pattern)

    # Build the upper evaluator
    ImpUpperEval = ImplicitUpperEvaluator()
    build_upper_evaluator!(ImpUpperEval, obj = f, constr = g, impfun = h,
                           nx = nx, np = np, nx = nx, ng = ng,
                           user_sparse = sparse_pattern)

    # Add nlp data blocks ("SHOULD" BE THE LAST THING TO DO)

    # Optimizes the model with load function
    MOI.optimize!(opt, CustomMod! = ImplicitMod, CustomModArgs = (ImpLowerEval,ImpUpperEval))

    return var,opt
end
