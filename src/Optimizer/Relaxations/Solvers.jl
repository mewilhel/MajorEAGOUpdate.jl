GetOptimizerType(x::RelaxationScheme) = x.OptimizerType

function GetOptimizer(x::Optimizer,y::Symbol)
    (y == :LP) && return x.LPOptimizer
    (y == :MILP) && return x.MILPOptimizer
    (y == :NLP) && return x.NLPOptimizer
    (y == :MINLP) && return x.MINLPOptimizer
    (y == :Custom1)  && return x.Custom1Optimizer
    (y == :Custom2)  && return x.Custom2Optimizer
    (y == :Custom3)  && return x.Custom3Optimizer
    error("Unrecognized solver type.")
end

# LINK TO ENVIRONMENT VARIABLES AND PARSE... NO DEPENDENCE ON SOLVERS DESIRABLE
function GetDefaultSolver(y::Symbol)
    if (y == :LP)
        error("No default LP solver set. Please use the SetDefaultLP(x::AbstractOptimizer)
               to set a default solver. Doing this once will store the default solver to
               your environment variables for future sessions.")
    elseif (y == :MILP)
        error("No default MILP solver set. Please use the SetDefaultMILP(x::AbstractOptimizer)
               to set a default solver. Doing this once will store the default solver to
               your environment variables for future sessions.")
    elseif (y == :NLP)
        error("No default NLP solver set. Please use the SetDefaultNLP(x::AbstractOptimizer)
               to set a default solver. Doing this once will store the default solver to
               your environment variables for future sessions.")
    elseif (y == :MINLP)
        error("No default MINLP solver set. Please use the SetDefaultMINLP(x::AbstractOptimizer)
               to set a default solver. Doing this once will store the default solver to
               your environment variables for future sessions.")
    else
        error("In order to use a default solver, your relaxation must specify :LP,
               :MILP, :NLP, or :MINLP. A custom solver should be set using the
               SetCustomOptimizer(x::AbstractOptimizer,i::Int) where i is in {1,2,3}.")
    end
end
