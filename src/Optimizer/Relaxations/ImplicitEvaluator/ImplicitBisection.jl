function ImplicitBisection(B::Optimizer,S::NodeBB)
  Pos = 0; Max = -Inf; TempMax = 0.0
  for i in 1:B.WorkingEvaluatorBlock.evaluator.np
    TempMax = (S.UpperVar[i] - S.LowerVar[i])/(B.VariableInfo[i].upper_bound - B.VariableInfo[i].lower_bound)
    if TempMax > Max
      Pos = i; Max = TempMax
    end
  end
  CutInterval = Interval(S.LowerVar[Pos],S.UpperVar[Pos])
  N1::Interval{Float64}, N2::Interval{Float64} = bisect(CutInterval)
  S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
  S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
  X1 = NodeBB(S.LowerVar, S.UpperVar, S.LowerBound, S.UpperBound, S.Depth + 1, -1, false)
  X2 = deepcopy(X1)
  X1.LowerVar[Pos] = N1.lo; X1.UpperVar[Pos] = N1.hi
  X2.LowerVar[Pos] = N2.lo; X2.UpperVar[Pos] = N2.hi
  return X1, X2
end
