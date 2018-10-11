DefaultCutCondition(x::Optimizer) = false
#CutContConditions(x::Optimizer) = (x.CutIterations < x.MaxCutIterations) ? false : CheckTightCut(x)
