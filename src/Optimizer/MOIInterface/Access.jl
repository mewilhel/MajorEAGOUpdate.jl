#=
MOI.canget(::EAGOOptimizer, ::MOI.NumberOfVariables) = true
MOI.canget(::EAGOOptimizer, ::MOI.ListOfVariableIndices) = true
MOI.canget(::EAGOOptimizer, ::MOI.SolverName) = true
MOI.canget(m::EAGOOptimizer, ::MOI.TerminationStatus) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveValue) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.ObjectiveBound) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.RelativeGap) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.SolveTime) = m.started_solve
MOI.canget(m::EAGOOptimizer, ::MOI.NodeCount) = m.started_solve
=#

MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [MOI.VariableIndex(i) for i in 1:length(m.VariableInfo)]
MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m.ObjectiveValue
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = length(m.VariableInfo)
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m.GlobalUpperBound
MOI.get(m::Optimizer, ::MOI.RelativeGap) = m.GlobalUpperBound
MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m.TerminationStatusCode
function MOI.get(m::Optimizer, ::MOI.SolveTime)
    IterCount = m.CurrentIterationCount
    if IterCount > 0
        return m.History.PreprocessTime[m.CurrentIterationCount] +
               m.History.PostprocessTime[m.CurrentIterationCount] +
               m.History.LowerTime[m.CurrentIterationCount] +
               m.History.UpperTime[m.CurrentIterationCount]
     else
         return 0.0
     end
 end
MOI.get(m::Optimizer, ::MOI.NodeCount) = length(m.MaximumNodeID)
