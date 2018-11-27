mutable struct VariableInfo
    is_integer::Bool
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
    start::Float64
end

# The default start value is zero.
VariableInfo() = VariableInfo(false,-Inf, false, Inf, false, false, 0.0)

LowerBound(x::VariableInfo) = x.lower_bound
UpperBound(x::VariableInfo) = x.upper_bound

@enum SolverFailed LowerSolverFailed UpperSolverFailed PreprocessingFailed PostProcessingFailed NoFailure
@enum OptimizerType LP MILP NLP MINLP

DummyFunction() = nothing

export Optimizer
mutable struct Optimizer <: MOI.AbstractOptimizer

    InputModel::Any
    IntegerVar::Vector{Int}
    VariableInfo::Vector{VariableInfo}
    LowerVariables::Vector{MOI.VariableIndex}
    UpperVariables::Vector{MOI.VariableIndex}
    VariableIndexLow::Vector{Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}}
    VariableIndexUpp::Vector{Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}}
    ObjConstrIndx::Vector{MOI.ConstraintIndex}
    InitialContinuousValues::IntervalBox            # Interval box constraints
    InitialIntegerValues::Vector{Int}               # Potential Integer Values

    NonlinearVariable::Dict{Int,Bool}
    PseudoCostLower::Vector{Float64}
    PseudoCostUpper::Vector{Float64}
    ProbCountLower::Vector{Float64}
    ProbCountUpper::Vector{Float64}

    VItoSto::Dict{Int,Int}
    StoToVI::Dict{Int,Int}
    ConstraintConvexity::Dict{MOI.ConstraintIndex,Bool}
    ConstraintLabel::Dict{Int,Symbol}

    ContinuousSolution::Vector{Float64}             # Stores a point in the IntervalSolutionBox
    IntegerSolution::Vector{Bool}                   # Stores the integer solution point
    Stack::Dict{Int,NodeBB}                       # Map of Node ID to NodeData

    NLPData::MOI.NLPBlockData
    WorkingEvaluatorBlock::MOI.NLPBlockData
    VariableNumber::Int
    ContinuousNumber::Int
    IntegerNumber::Int
    ConstraintNumber::Int
    LinearNumber::Int
    QuadraticNumber::Int

    CurrentLowerInfo::LowerInfo                       # Problem solution info for lower bounding program
    CurrentUpperInfo::UpperInfo                       # Problem solution info for upper bounding program
    CurrentPreprocessInfo::PreprocessInfo             # Problem solution info for preprocessing step
    CurrentPostprocessInfo::PostprocessInfo           # Problem solution info for postprocessing step

    Objective::Union{MOI.SingleVariable,MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64},Nothing}
    ObjectiveConvexity::Bool

    LinearLEQConstraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64},Int}}
    LinearGEQConstraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64},Int}}
    LinearEQConstraints::Vector{Tuple{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64},Int}}
    QuadraticLEQConstraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64},Int}}
    QuadraticGEQConstraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64},Int}}
    QuadraticEQConstraints::Vector{Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64},Int}}
    QuadraticConvexity::Vector{Bool}

    UniQuadraticLEQConstraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    UniQuadraticGEQConstraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    UniQuadraticEQConstraints::Vector{Tuple{Float64,Float64,Float64,Int}}
    BiQuadraticLEQConstraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}
    BiQuadraticGEQConstraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}
    BiQuadraticEQConstraints::Vector{Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}}

    LPOptimizer::MOI.AbstractOptimizer
    MILPOptimizer::MOI.AbstractOptimizer
    NLPOptimizer::MOI.AbstractOptimizer
    MINLPOptimizer::MOI.AbstractOptimizer

    InitialRelaxedOptimizer::MOI.AbstractOptimizer
    WorkingRelaxedOptimizer::MOI.AbstractOptimizer
    InitialUpperOptimizer::MOI.AbstractOptimizer
    WorkingUpperOptimizer::MOI.AbstractOptimizer

    Relaxation::RelaxationScheme

    # Stores functions for implementing Branch and Bound
    LowerProblem::Function                                                 # Stores lower problem function
    UpperProblem::Function                                                 # Stores upper problem function
    Preprocess::Function                                                   # Preprocessing function
    Postprocess::Function                                                  # Post processing function
    RepeatCheck::Function                                                  # Repeation check
    ConvergenceCheck::Function                                             # convergence criterion
    TerminationCheck::Function                                             # Stores termination check function
    NodeStorage::Function                                                  # Stores branching function
    NodeSelection::Function                                                # Stores node selection function
    BisectionFunction::Function                                            #
    CutCondition::Function                                                 #
    AddCut::Function                                                       #

    GlobalLowerBound::Float64                                              # Global Lower Bound
    GlobalUpperBound::Float64                                              # Global Upper Bound
    MaximumNodeID::Int
    History::NodeHistory
    CurrentIterationCount::Int
    CurrentNodeCount::Int

    # Storage for output
    SolutionValue::Float64                                                 # Value of the solution
    FirstFound::Bool                                                       #
    FeasibleSolutionFnd::Bool                                              #
    FirstSolutionNode::Int                                                 #
    LastGap::Float64                                                       #
    OptimizationSense::MOI.OptimizationSense                               #
    ObjectiveValue::Float64
    TerminationStatusCode::MOI.TerminationStatusCode                       #
    ResultStatusCode::MOI.ResultStatusCode
    StartedSolve::Bool                                                     #
    FailedSolver::SolverFailed                                             # Stores branching function

    # Output specification fields
    Verbosity::Int                                                         # Stores output selection
    WarmStart::Bool                                                        # Boolean
    OutputInterations::Int                                                 # Number of iterations to skip between printing iteration summary
    HeaderInterations::Int                                                 # Number of iterations to skip between printing heade
    DigitsDisplayed::Int                                                   # digits displayed before decimal
    ReturnHistory::Bool                                                    # returns LBD, UBD array and time vector
    FlagSubSolverErrors::Bool                                              # If a subsolver has a problem termination code then stop the algorithm
                                                                           # and record it
    # Termination Limits
    IterationLimit::Int                                                    # Maximum number of iterations
    NodeLimit::Int                                                         # Maximum number of nodes to store in memory
    AbsoluteTolerance::Float64                                             # Absolute tolerance for BnB
    RelativeTolerance::Float64                                             # Relative tolerance for BnB
    ExhaustiveSearch::Bool                                                 # Exhaustive search: find all solns or find first

    # Optimality-Based Bound Tightening (OBBT) Options
    OBBTVars::Vector{MOI.VariableIndex}
    OBBTDepth::Int
    OBBTRepts::Int
    OBBTAggrOn::Bool
    OBBTAggrMaxIteration::Int
    OBBTAggrMinDimLimit::Int
    OBBTTolerance::Float64
    OBBTWorkingLowerIndx::Vector{MOI.VariableIndex}
    OBBTWorkingUpperIndx::Vector{MOI.VariableIndex}
    OBBTInitialLowerIndx::Vector{MOI.VariableIndex}
    OBBTInitialUpperIndx::Vector{MOI.VariableIndex}
    OBBTActiveCurrent::Bool

    # Duality-Based Bound Tightening (DBBT) Options
    DBBTDepth::Int
    DBBTTolerance::Float64

    # Feasibility-Based Bound Tightening Options
    CPWalkDepth::Int
    CPWalkRepts::Int
    EvalWalkRepts::Int
    EvalReverse::Bool

    # Options for Poor Man's LP
    PoorManLPDepth::Int
    PoorManLPRepts::Int

    # Options for Quadratic Bounding Tightening
    UniQuadDepth::Int
    UniQuadRepts::Int
    BiQuadDepth::Int
    BiQuadRepts::Int

    # Options for Repetition (If DBBT Performed Well)
    MaximumRepetitions::Int
    RepetitionVolumeTolerance::Float64

    # Cutting Plane Options
    CutIterations::Int

    # Status flags
    FirstRelaxed::Bool

    # Debug
    Debug1::Any
    Debug2::Any

    function Optimizer(;options...)

        m = new()

        default_opt_dict = Dict{Symbol,Any}()

        # set fallback for potentially user defined functions
        for i in (:LowerProblem, :UpperProblem, :Preprocess, :Postprocess, :RepeatCheck,
                  :ConvergenceCheck, :TerminationCheck, :NodeStorage, :NodeSelection,
                  :BisectionFunction, :CutCondition, :AddCut)
                  default_opt_dict[i] = DummyFunction
        end

        # set fallback for optimizers
        for i in (:LPOptimizer, :MILPOptimizer, :NLPOptimizer, :MINLPOptimizer,
                  :InitialRelaxedOptimizer, :WorkingRelaxedOptimizer,
                  :InitialUpperOptimizer, :WorkingUpperOptimizer)
            default_opt_dict[i] = DummyOptimizer()
        end

        default_opt_dict[:Relaxation] = DefaultRelaxationScheme()
        default_opt_dict[:GlobalLowerBound] = -Inf
        default_opt_dict[:GlobalUpperBound] = Inf

        # Output specification fields
        #m.Verbosity = 0
        default_opt_dict[:Verbosity] = 3
        default_opt_dict[:WarmStart] = false
        default_opt_dict[:OutputInterations] = 1
        default_opt_dict[:HeaderInterations] = 100
        default_opt_dict[:DigitsDisplayed] = 3
        default_opt_dict[:ReturnHistory] = false
        default_opt_dict[:FlagSubSolverErrors] = true

        # Duality-based bound tightening parameters
        default_opt_dict[:DBBTDepth] = Int(1E6)
        default_opt_dict[:DBBTTolerance] = 1E-8

        # Optimality-based bound tightening parameters
        default_opt_dict[:OBBTDepth] = 3
        default_opt_dict[:OBBTAggrOn] = false
        default_opt_dict[:OBBTAggrMaxIteration] = 2
        default_opt_dict[:OBBTAggrMinDimLimit] = 2
        default_opt_dict[:OBBTTolerance] = 1E-5

        # Feasibility-Based Bound Tightening Options
        default_opt_dict[:CPWalkDepth] = 10
        default_opt_dict[:CPWalkRepts] = 10
        default_opt_dict[:EvalWalkRepts] = 1
        default_opt_dict[:EvalReverse] = false

        # Options for Repetition (If DBBT Performed Well)
        default_opt_dict[:MaximumRepetitions] = 1
        default_opt_dict[:RepetitionVolumeTolerance] = 0.0

        # Poor Man's LP reptiations
        default_opt_dict[:PoorManLPRepts] = 1

        # Termination Limits
        default_opt_dict[:IterationLimit] = Int(1E6)
        default_opt_dict[:NodeLimit] = Int(1E6)
        default_opt_dict[:AbsoluteTolerance] = 1E-4
        default_opt_dict[:RelativeTolerance] = 1E-4
        default_opt_dict[:ExhaustiveSearch] = false
        default_opt_dict[:FirstRelaxed] = false

        for i in keys(default_opt_dict)
            if (haskey(options,i))
                setfield!(m, i, options[i])
            else
                setfield!(m, i, default_opt_dict[i])
            end
        end

        m.Debug1 = []
        m.Debug2 = []
        m.InputModel = 0
        m.IntegerVar = Int[]
        m.VariableInfo = VariableInfo[]
        m.LowerVariables = MOI.VariableIndex[]
        m.UpperVariables = MOI.VariableIndex[]
        m.VariableIndexLow = Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}[]
        m.VariableIndexUpp = Tuple{MOI.ConstraintIndex,MOI.ConstraintIndex,Int}[]
        m.ObjConstrIndx = MOI.ConstraintIndex[]
        m.InitialContinuousValues = IntervalBox(Interval(0.0))      # Interval box constraints
        m.InitialIntegerValues = Vector{Int}[]                      # Potential Integer Values
        m.NLPData = empty_nlp_data()
        m.WorkingEvaluatorBlock = empty_nlp_data()

        m.ConstraintConvexity = Dict{MOI.ConstraintIndex,Bool}()
        m.VItoSto = Dict{Int,Int}()
        m.StoToVI = Dict{Int,Int}()
        m.NonlinearVariable = Dict{Int,Bool}()
        m.PseudoCostLower = Float64[]
        m.PseudoCostUpper = Float64[]
        m.ProbCountLower = Float64[]
        m.ProbCountUpper = Float64[]

        m.ContinuousSolution = Float64[]                            # Stores a point in the IntervalSolutionBox
        m.IntegerSolution = Bool[]                                  # Stores the integer solution point
        m.Stack = Dict{Int,NodeBB}()                              # Map of Node ID to NodeData
        m.VariableNumber = 0
        m.ContinuousNumber = 0
        m.IntegerNumber = 0
        m.ConstraintNumber = 0

        m.LinearLEQConstraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}, Int}[]
        m.LinearGEQConstraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, Int}[]
        m.LinearEQConstraints = Tuple{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}, Int}[]
        m.QuadraticLEQConstraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.LessThan{Float64}, Int}[]
        m.QuadraticGEQConstraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.GreaterThan{Float64}, Int}[]
        m.QuadraticEQConstraints = Tuple{MOI.ScalarQuadraticFunction{Float64}, MOI.EqualTo{Float64}, Int}[]

        m.UniQuadraticLEQConstraints = Tuple{Float64,Float64,Float64,Int}[]
        m.UniQuadraticGEQConstraints = Tuple{Float64,Float64,Float64,Int}[]
        m.UniQuadraticEQConstraints = Tuple{Float64,Float64,Float64,Int}[]
        m.BiQuadraticLEQConstraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]
        m.BiQuadraticGEQConstraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]
        m.BiQuadraticEQConstraints = Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int,Int}[]

        m.CurrentLowerInfo = LowerInfo()
        m.CurrentUpperInfo = UpperInfo()
        m.CurrentPreprocessInfo = PreprocessInfo()
        m.CurrentPostprocessInfo = PostprocessInfo()

        m.Objective = nothing
        m.ObjectiveConvexity = false

        # Historical Information
        m.History = NodeHistory()
        m.CurrentIterationCount = 0
        m.CurrentNodeCount = 0
        m.MaximumNodeID = 0

        # Output for Solution Storage
        m.SolutionValue = -Inf
        m.FirstFound = false
        m.FeasibleSolutionFnd = false
        m.FirstSolutionNode = -1
        m.LastGap = -Inf
        m.OptimizationSense = MOI.FeasibilitySense
        m.TerminationStatusCode = MOI.OtherError #MOI.OptimizeNotCalled
        m.ResultStatusCode = MOI.OtherResultStatus
        m.StartedSolve = false
        m.FailedSolver = NoFailure

        # Optimality-Based Bound Tightening (OBBT) Storage
        m.OBBTVars = MOI.VariableIndex[]
        m.OBBTWorkingLowerIndx = MOI.VariableIndex[]
        m.OBBTWorkingUpperIndx = MOI.VariableIndex[]
        m.OBBTInitialLowerIndx = MOI.VariableIndex[]
        m.OBBTInitialUpperIndx = MOI.VariableIndex[]

        return m
    end

end
