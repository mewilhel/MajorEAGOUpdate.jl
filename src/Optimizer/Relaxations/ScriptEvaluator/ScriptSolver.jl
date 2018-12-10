# Modifies functions post initial relaxation to use appropriate nlp evaluator
function ScriptMod(opt::Optimizer,args)
    ImpLowerEval = args[1]; ImpUpperEval = args[2]; has_obj = args[3]; c_bnds = args[4]
    lower_eval_block = MOI.NLPBlockData(c_bnds, ImpLowerEval, has_obj)
    opt.RelaxFunction = ImplicitRelaxModel!
    opt.WorkingEvaluatorBlock = deepcopy(lower_eval_block)
    # load lower nlp data block
    if MOI.supports(opt.InitialRelaxedOptimizer, MOI.NLPBlock())
        if ~isempty(opt.NonlinearVariable)
            opt.InitialRelaxedOptimizer.nlp_data = deepcopy(lower_eval_block)
            opt.WorkingRelaxedOptimizer.nlp_data = deepcopy(lower_eval_block)
        end
    end
    # load upper nlp data block
    upper_eval_block = MOI.NLPBlockData(c_bnds, ImpUpperEval, has_obj)
    if MOI.supports(opt.InitialRelaxedOptimizer, MOI.NLPBlock())
        if ~isempty(opt.NonlinearVariable)
            opt.InitialRelaxedOptimizer.nlp_data = deepcopy(upper_eval_block)
            opt.WorkingRelaxedOptimizer.nlp_data = deepcopy(upper_eval_block)
        end
    end
end

function BuildScriptEval(f,g,h)
end

function SolveScript(f, h, xl, xu, pl, pu, opt, hj, g)

    # get dimensions
    @assert length(pl) == length(pu)
    @assert length(xl) == length(xu)
    np = length(pl); nx = length(xl);
    ng = (g == nothing) ? 0 : length(g(xl,pl))

    # sets most routines to default (expect bisection)
    SetToDefault!(opt)

    # add variables to lower, upper, and EAGO models
    var_EAGO = MOI.add_variables(opt, np+nx)

    for j in 1:nx
        MOI.add_constraint(opt, var_EAGO[j], MOI.GreaterThan(xl[j]))
        MOI.add_constraint(opt, var_EAGO[j], MOI.LessThan(xu[j]))
    end

    for i in 1:np
        MOI.add_constraint(opt, var_EAGO[i+nx], MOI.GreaterThan(pl[i]))
        MOI.add_constraint(opt, var_EAGO[i+nx], MOI.LessThan(pu[i]))
    end

    # Build the lower implicit evaluator

    # Build the upper evaluator
    ImpUpperEval = ImplicitUpperEvaluator()
    build_upper_evaluator!(ImpUpperEval, h, np, nx, obj = f, constr = g, ng = ng) #, user_sparse = user_sparsity)

    # Add nlp data blocks ("SHOULD" BE THE LAST THING TO DO)
    has_obj = (f != nothing)
    bnd_pair = MOI.NLPBoundsPair(-Inf,0.0)
    nlp_bnds = [bnd_pair for i=1:ng]

    # Optimizes the model with load function
    MOI.optimize!(opt, CustomMod! = ScriptMod, CustomModArgs = (ImpLowerEval,ImpUpperEval,has_obj,nlp_bnds))

    return var_EAGO,opt
end
#=
  g = DummyFunction, hj = DummyFunction,
  user_sparsity = Tuple{Int64,Int64}[])
=#
