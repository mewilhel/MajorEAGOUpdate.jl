# Descibe on intervalbox versus interval vector

function ContinuousRelativeBisect(B::Optimizer,S::NodeBB)
  Pos = indmax(S.UpperVars - S.LowerVars)
  N1::Interval{Float64},N2::Interval{Float64} = bisect(Interval(S.LowerVars[Pos],S.UpperVars[Pos]))
  S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
  S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
  X1 = copy(S); X1.LowerVars[Pos] = Nl.lo; X1.UpperVars[Pos] = Nl.hi; X1.Depth += 1
  X2 = copy(S); X2.LowerVars[Pos] = N2.lo; X2.UpperVars[Pos] = N2.hi; X2.Depth += 1
  return X1,X2
end

#
function PseudoCostBranch(B::Optimizer,S::NodeBB)
    indx, psval = -1,-Inf
    for var in B.NonlinearVariable
      LowerCount = B.ProbCountLower[var]
      CostLower =(LowerCount > 0.0) ? (B.PseudoCostLower[var]/LowerCount) : 0.0
      UpperCount = B.ProbCountLower[var]
      CostUpper =(UpperCount > 0.0) ? (B.PseudoCostUpper[var]/UpperCount) : 0.0
      Score = max(CostLower,B.PseudoTol)*min(CostUpper,B.PseudoTol)
      if (Score > psval)
          psval = Score
          indx = var
      end
    end
    B.ProbCountLower[indx] += 1; B.ProbCountLower[var]
    return index, psval
end

# hybrid reliability rule branching on integer variables,
function PseudoCostBisect(B::Optimizer,S::NodeBB)
  # if fractional integer terms --> branch on them (hybrid reliability rule)
  FractionalInteger = is_integer_feasible(B)
  Pos,Value = PseudoCostBranch(B,S)
  if FractionalInteger
      S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
      S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
      X1 = copy(S); X1.LowerVars[Pos] = 0.0; X1.UpperVars[Pos] = 0.0; X1.Depth += 1; X1.LastBranch = Pos; X1.DirBranch = false
      X2 = copy(S); X2.LowerVars[Pos] = 1.0; X2.UpperVars[Pos] = 1.0; X2.Depth += 1; X2.LastBranch = Pos; X2.DirBranch = true
  else
      N1::Interval{Float64},N2::Interval{Float64} = bisect(Interval(S.LowerVars[Pos],S.UpperVars[Pos]),Value)
      S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
      S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
      X1 = copy(S); X1.LowerVars[Pos] = Nl.lo; X1.UpperVars[Pos] = Nl.hi; X1.Depth += 1; X1.LastBranch = Pos; X1.DirBranch = false
      X2 = copy(S); X2.LowerVars[Pos] = N2.lo; X2.UpperVars[Pos] = N2.hi; X2.Depth += 1; X2.LastBranch = Pos; X2.DirBranch = true
  end
end
