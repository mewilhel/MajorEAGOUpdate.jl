function TrivFunction(x) end

function MOI.optimize!(m::Optimizer; CustomMod! = TrivFunction, CustomModArgs = (1,))

    ########### Reformulate DAG using auxilliary variables ###########
    #LoadDAG!(m); LabelDAG!(m)
    #NewVariableSize,NewVariableIndex = ProcessDAG!(m)
    NewVariableSize = length(m.VariableInfo)
    m.ContinuousNumber = NewVariableSize
    m.VariableNumber = NewVariableSize

    ########### Set Correct Size for Problem Storage #########
    m.CurrentLowerInfo.Solution = Float64[0.0 for i=1:NewVariableSize]
    m.CurrentLowerInfo.LowerVarDual = Float64[0.0 for i=1:NewVariableSize]
    m.CurrentLowerInfo.UpperVarDual = Float64[0.0 for i=1:NewVariableSize]
    m.CurrentUpperInfo.Solution = Float64[0.0 for i=1:NewVariableSize]

    # loads variables into the working model
    m.VItoSto = Dict{Int,Int}()
    for i=1:NewVariableSize
        m.VItoSto[i] = i
    end
    m.StoToVI = ReverseDict(m.VItoSto)

    ###### OBBT Setup #####
    # Sets terms that OBBT will be performed on
    for i=1:NewVariableSize
        if m.NonlinearVariable[i]
            push!(m.OBBTVars,MOI.VariableIndex(i))
        end
    end

    # Get various other sizes
    num_nlp_constraints = length(m.NLPData.constraint_bounds)
    m.ContinuousSolution = zeros(Float64,NewVariableSize)

    # Sets any unset functions to default values
    SetToDefault!(m)

    # Copies variables and bounds to lower subproblems
    #PushLowerVariables!(m)
    m.UpperVariables = MOI.add_variables(m.InitialRelaxedOptimizer, m.VariableNumber)

    # Create initial node and add it to the stack
    CreateInitialNode!(m)

    # Build the JuMP NLP evaluator
    evaluator = m.NLPData.evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad]
    #has_hessian && push!(init_feat, :Hess)
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator,init_feat)

    # Sets up relaxations terms that don't vary during iterations (mainly linear)
    m.LowerVariables = MOI.VariableIndex.(1:m.VariableNumber)

    # Check convexity of objective function and quadratic constraints

    # Creates initial nlp evaluator
    m.WorkingEvaluatorBlock = m.NLPData
    if (typeof(m.NLPData.evaluator) != MajorEAGOUpdate.EmptyNLPEvaluator)
        Built_Evaluator = Build_NLP_Evaluator(MC{m.VariableNumber},m.NLPData.evaluator,m)
        (m.OptimizationSense == MOI.MaxSense) && MinusObjective!(Built_Evaluator)
        m.WorkingEvaluatorBlock = MOI.NLPBlockData(m.NLPData.constraint_bounds, Built_Evaluator, m.NLPData.has_objective)
    end

    # Relax initial model terms
    RelaxModel!(m, m.InitialRelaxedOptimizer, m.Stack[1], m.Relaxation, load = true)

    # Runs a customized function if one is provided
    m.CustomModFlag = (CustomMod! != TrivFunction)
    if m.CustomModFlag
        CustomMod!(m,CustomModArgs)
    end

    # if optimizer type is supplied for upper, build factory
    if (~m.UseUpperFactory)
        MOI.add_variables(m.InitialUpperOptimizer, m.VariableNumber)
        m.UpperVariables = MOI.add_variables(m.WorkingUpperOptimizer, m.VariableNumber)
        SetLocalNLP!(m)
    end

    # Runs the branch and bound routine
    SolveNLP!(m)
end
