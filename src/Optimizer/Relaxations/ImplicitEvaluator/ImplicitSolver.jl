function SolveImplicit(f::Function, g::Function, h::Function,
                       xl::Vector{Float64}, xu::Vector{Float64},
                       pl::Vector{Float64}, pu::Vector{Float64},
                       opt::Optimizer; user_sparsity::Vector{Tuple{Int64,Int64}} = sparse_pattern)

    #
    @assert length(pl) == length(pu)
    @assert length(xl) == length(xu)
    np = length(pl); nx = length(xl)

    # Sets most routines to default
    SetToDefault!(opt)
    opt.BisectionFunction = ImplicitBisection

    # add variables to upper and lower models
    MOI.add_variables(opt.InitialLowerOptimizer, np+nx)
    MOI.add_variables(opt.InitialUpperOptimizer, np+nx)

    # Build the lower evaluator
    ImpLowerEval = ImplicitLowerEvaluator()
    build_lower_evaluator!(ImpLowerEval, obj = f, constr = g, impfun = h,
                           nx = nx, np = np, user_sparse = sparse_pattern)

    # Build the upper evaluator
    ImpUpperEval = ImplicitUpperEvaluator()
    ImpUpperEval.np = np

    # Creates the initial node
    opt.Stack[1] = NodeBB()
    opt.Stack[1].LowerVar = vcat(pl,xl)
    opt.Stack[1].UpperVar = vcat(pu,xu)
    opt.MaximumNodeID += 1

    # creates the Implicit Evaluator
    # loads the Implicit Evaluator into the optimizer
    # solves the problem
    # outputs the results
end
