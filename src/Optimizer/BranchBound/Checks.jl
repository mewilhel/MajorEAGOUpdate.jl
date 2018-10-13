```
    DefaultRepeatCheck

ADD DOC STRING
```
function DefaultRepeatCheck(x::Optimizer,y::NodeBB)
  return false
end

```
    DefaultTerminationCheck

ADD DOC STRING
```
function DefaultTerminationCheck(x::Optimizer)

    L = x.CurrentLowerInfo.Value
    U = x.CurrentUpperInfo.Value

    t1 = length(x.Stack) > 0                                   # stack is non-empty
    t2 = x.CurrentIterationCount < x.IterationLimit            # maximum iterations not exceeded
    t3 = length(x.Stack) < x.NodeLimit                         # maximum node number not exceeded
    t4 = (U - L) > x.AbsoluteTolerance                         # absolute tolerance satisfied
    t5 = (abs(U - L)/(min(abs(L),abs(U))) > x.RelativeTolerance) || ~(L > -Inf)   # relative tolerance satisfied

    if t1 & t2 & t3 & t4 & t5
      return true
    else
      if (x.Verbosity == 2 || x.Verbosity == 3)
        if ~t1
          if (x.FirstSolutionNode > 0)
            println("Empty Stack: Exhaustive Search Finished")
          else
            println("Empty Stack: Infeasible")
          end
        elseif ~t3
          println("Node Limit Exceeded")
        elseif ~t2
          println("Maximum Iteration Exceeded")
        else
          println("Convergence Tolerance Reached")
        end
      end
      return false
    end
end

```
    DefaultConvergenceCheck

ADD DOC STRING
```
function DefaultConvergenceCheck(x::Optimizer)
  L = x.CurrentLowerInfo.Value
  U = x.CurrentUpperInfo.Value
  t1 = (U - L) <= x.AbsoluteTolerance
  t2 = (abs(U - L)/(min(abs(L),abs(U))) <= x.RelativeTolerance)
  return t1 || t2
end
