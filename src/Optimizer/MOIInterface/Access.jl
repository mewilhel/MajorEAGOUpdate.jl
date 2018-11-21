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
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m.ResultStatusCode

function MOI.get(m::Optimizer, ::MOI.SolveTime)
    if m.CurrentIterationCount > 0
        return m.History.PreprocessTime[m.CurrentIterationCount-1] +
               m.History.PostprocessTime[m.CurrentIterationCount-1] +
               m.History.LowerTime[m.CurrentIterationCount-1] +
               m.History.UpperTime[m.CurrentIterationCount-1]
     else
         return 0.0
     end
 end
MOI.get(m::Optimizer, ::MOI.NodeCount) = length(m.MaximumNodeID)

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return model.ContinuousSolution[vi.value]
end
