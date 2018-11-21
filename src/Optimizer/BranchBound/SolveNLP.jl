# () Iteration Printing

```
    SolveNLP!

Solves the branch and bound problem with the input model and solver object.
```
function SolveNLP!(x::Optimizer)

  # initializes Flags
  TT = stdout
  UBDSaveFlag = false
  UBDTempSaveFlag = false
  PostSaveFlag = false
  if (~x.WarmStart)
    x.CurrentIterationCount = 1
    x.CurrentNodeCount = 0
  end
  NumberOfCuts = 0

  # terminates when max nodes or iteration is reach, or when node stack is empty
  iterationcountinternal = 0
  while (x.TerminationCheck(x))
    iterationcountinternal += 1
    #println("iterationcountinternal: $iterationcountinternal")
    # Fathom nodes with lower bound greater than global upper bound
    x.GlobalLowerBound = FindLowerBound(x)
    x.History.LowerBound[x.CurrentIterationCount] = x.GlobalLowerBound
    # Sets conditional save flags
    UBDSaveFlag = false
    UBDTempSaveFlag = false
    PostSaveFlag = false
    x.CutIterations = 0

     # Selects node, deletes it from stack, prints based on verbosity
    CurrentKey,CurrentNode = x.NodeSelection(x)
    (x.Verbosity == 3) && PrintNode!(CurrentKey,CurrentNode) # Prints node in full verbosity mode

    # Solves preprocessing/LBD/UBD/postprocessing once to get timing right
    x.CurrentPreprocessInfo.Feasibility = true
    x.CurrentPostprocessInfo.Feasibility = true
    if (x.CurrentIterationCount == 1)
      #println("start preprocess")
      OldLowerInfo = deepcopy(x.CurrentLowerInfo)
      OldUpperInfo = deepcopy(x.CurrentUpperInfo)
      OldPreprocessInfo = deepcopy(x.CurrentPreprocessInfo)
      OldPostprocessInfo = deepcopy(x.CurrentPostprocessInfo)
      tempNode = copy(CurrentNode)
      x.Preprocess(x,tempNode)
      #println("ran initial preprocess")
      #x.LowerProblem(x,tempNode)
      #println("ran initial lower problem")
      x.UpperProblem(x,tempNode)
      #println("ran initial upper problem")
      x.Postprocess(x,tempNode)
      #println("ran initial postprocess")
      x.CurrentLowerInfo = OldLowerInfo
      x.CurrentUpperInfo = OldUpperInfo
      x.CurrentPreprocessInfo = OldPreprocessInfo
      x.CurrentPostprocessInfo = OldPostprocessInfo
    end
    x.CurrentPreprocessInfo.Feasibility = true
    x.CurrentPostprocessInfo.Feasibility = true

    # Performs prepocessing and times
    redirect_stdout()
    PreprocessTime = @elapsed x.Preprocess(x,CurrentNode)
    PreprocessTime = x.Preprocess(x,CurrentNode)
    x.History.PreprocessTime[x.CurrentIterationCount] = x.History.PreprocessTime[x.CurrentIterationCount-1]+PreprocessTime
    redirect_stdout(TT)

    x.CurrentUpperInfo.Feasibility = true
    if (x.CurrentPreprocessInfo.Feasibility)
      # solves & times lower bounding problem

      x.LowerProblem(x,CurrentNode)
      LowerProblemTime = 0.0
      x.History.LowerBound[x.CurrentIterationCount] = x.History.LowerBound[x.CurrentIterationCount-1]+LowerProblemTime
      x.History.LowerTime[x.CurrentIterationCount] = x.History.LowerTime[x.CurrentIterationCount-1]+LowerProblemTime
      x.History.LowerCount += 1
      PrintResults!(x,true)
      while x.CutCondition(x)
        x.AddCut!(x); x.CutIterations += 1
        redirect_stdout()
        #LowerProblemTime = @elapsed x.LowerProblem(x,CurrentNode)
        x.LowerProblem(x,CurrentNode)
        LowerProblemTime = 0.0
        x.History.LowerBound[x.CurrentIterationCount] = x.History.LowerBound[x.CurrentIterationCount]+LowerProblemTime
        x.History.LowerTime[x.CurrentIterationCount] = x.History.LowerTime[x.CurrentIterationCount-1]+LowerProblemTime
        x.History.LowerCount += 1
        redirect_stdout(TT)
        PrintResults!(x,true)
      end
      x.History.CutCount[x.CurrentIterationCount] = x.CutIterations
      println("post cut, lower feas: $(x.CurrentLowerInfo.Feasibility)")
      # checks for infeasibility stores solution
      if (x.CurrentLowerInfo.Feasibility)
        println("lower is feas")
        if (~x.ConvergenceCheck(x))
          println("didn't converge")
          println("CurrentNode: $(CurrentNode)")
          # Solves upper bounding problem
          #UpperProblemTime = @elapsed x.UpperProblem(x,CurrentNode)
          UpperProblemTime = 0.0
          x.UpperProblem(x,CurrentNode)
          x.History.UpperTime[x.CurrentIterationCount] = x.History.UpperTime[x.CurrentIterationCount-1]+UpperProblemTime
          x.History.UpperCount += 1
          UBDTempSaveFlag = true
          PrintResults!(x,false)

          println("Current Upper Value: $(x.CurrentUpperInfo.Value)")
          println("Current Upper Feasiblity: $(x.CurrentUpperInfo.Feasibility)")
          # Stores information if better feasible upper bound is formed
          if (x.CurrentUpperInfo.Feasibility)
            if (x.CurrentUpperInfo.Value < x.GlobalUpperBound)
              x.FeasibleSolutionFnd = true
              x.FirstSolutionNode = x.MaximumNodeID
              x.SolutionValue = x.CurrentUpperInfo.Value
              x.ContinuousSolution = x.CurrentUpperInfo.Solution
              x.History.UpperBound[x.CurrentIterationCount] = x.SolutionValue
              UBDSaveFlag = true
            end
          end

          # Performs and times post processing
          redirect_stdout()
          PostprocessTime = @elapsed x.Postprocess(x,CurrentNode)
          x.History.PostprocessTime[x.CurrentIterationCount] = x.History.PostprocessTime[x.CurrentIterationCount-1]+PostprocessTime
          PostSaveFlag = true
          redirect_stdout(TT)

          # Checks to see if the node
          if (x.CurrentPostprocessInfo.Feasibility)
            if x.RepeatCheck(x,CurrentNode)
              SingleStorage!(x,CurrentNode)
            else
              Y1,Y2 = x.BisectionFunction(x,CurrentNode)
              x.NodeStorage(x,Y1,Y2)
            end
          end
        elseif (~x.ExhaustiveSearch && x.FeasibleSolutionFnd)
        end
      else
        x.History.UpperBound[x.CurrentIterationCount] = x.History.UpperBound[x.CurrentIterationCount-1]
      end
      Fathom!(x)
    else
      x.CurrentLowerInfo.Value = -Inf
      x.CurrentLowerInfo.Feasibility = false
      x.CurrentUpperInfo.Feasibility = false
    end
    (~UBDSaveFlag) && (x.History.UpperBound[x.CurrentIterationCount] = x.GlobalUpperBound)
    (~UBDTempSaveFlag) && (x.History.UpperTime[x.CurrentIterationCount] = x.History.UpperTime[x.CurrentIterationCount-1])
    (~PostSaveFlag) && (x.History.PostprocessTime[x.CurrentIterationCount] = x.History.PostprocessTime[x.CurrentIterationCount-1])
    x.History.Count[x.CurrentIterationCount] = x.CurrentNodeCount
    PrintIteration!(x)
    x.CurrentIterationCount += 1
  end

  x.SolutionValue = x.GlobalUpperBound           # Sets the solution value found
  PrintSolution!(x)                              # Prints the solution
end
