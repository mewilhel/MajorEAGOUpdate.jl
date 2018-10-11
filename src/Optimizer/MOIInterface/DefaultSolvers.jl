

function SetDefaultLPSolver!(x::MOI.AbstractOptimizer)
    ENV["EAGO_LP_Solver"] = x
end
SetDefaultMILPSolver!(x::MOI.AbstractOptimizer) && (ENV["EAGO_MILP_Solver"] = x)
SetDefaultNLPSolver!(x::MOI.AbstractOptimizer) && (ENV["EAGO_NLP_Solver"] = x)
SetDefaultMINLPSolver!(x::MOI.AbstractOptimizer) && (ENV["EAGO_MINLP_Solver"] = x)

function GetOptimizer(x::Optimizer)

end

function GetDefaultSolver(x::Optimizer)
end
