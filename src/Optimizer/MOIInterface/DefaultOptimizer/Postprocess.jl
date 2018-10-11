function EAGODefault_PostProcess!(x::Optimizer,y::NodeBB)
    Variable_DBBT!(y,x.CurrentLowerInfo.LowerVarDual,
                     x.CurrentLowerInfo.UpperVarDual,
                     x.CurrentLowerInfo.Value,
                     x.GlobalUpperBound)
    x.CurrentPostprocessInfo.Feasibility = true
end
