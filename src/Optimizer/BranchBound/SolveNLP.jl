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
    x.CurrentIterationCount = 0
    x.CurrentNodeCount = 0
  end
  NumberOfCuts = 0

  # terminates when max nodes or iteration is reach, or when node stack is empty
  while (x.TerminationCheck(x))

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
    (x.Verbosity == :Full) && PrintNode!(CurrentKey,CurrentNode) # Prints node in full verbosity mode

    # Solves preprocessing/LBD/UBD/postprocessing once to get timing right
    x.CurrentPreprocessInfo.Feasibility = true
    x.CurrentPostprocessInfo.Feasibility = true
    if (x.CurrentIterationCount == 0)
      tempNode = copy(CurrentNode)
      x.Preprocess(x,tempNode)
      x.LowerProblem(x,tempNode)
      x.UpperProblem(x,tempNode)
      x.Postprocess(x,tempNode)
    end
    x.CurrentPreprocessInfo.Feasibility = true
    x.CurrentPostprocessInfo.Feasibility = true

    # Performs prepocessing and times
    #redirect_stdout()
    #PreprocessTime = @elapsed x.Preprocess(x,CurrentNode)
    PreprocessTime = x.Preprocess(x,CurrentNode)
    #x.History.PreprocessTime[x.CurrentIterationCount] = x.History.PreprocessTime[x.CurrentIterationCount-1]+PreprocessTime
    #redirect_stdout(TT)

    x.CurrentUpperInfo.Feasibility = true
    if (x.CurrentPreprocessInfo.Feasibility)
      # solves & times lower bounding problem

      #redirect_stdout()
      #LowerProblemTime = @elapsed x.LowerProblem(x,CurrentNode)
      LowerProblemTime = x.LowerProblem(x,CurrentNode)
      #x.History.LowerBound[x.CurrentIterationCount] = x.History.LowerBound[x.CurrentIterationCount-1]+LowerProblemTime
      x.History.LowerCount += 1
      #redirect_stdout(TT)
      PrintResults!(x,true)

      while x.CutCondition(x)
        x.AddCut!(x); x.CutIterations += 1
        #redirect_stdout()
        #LowerProblemTime = @elapsed x.LowerProblem(x,CurrentNode)
        LowerProblemTime = x.LowerProblem(x,CurrentNode)
        #x.History.LowerBound[x.CurrentIterationCount] = x.History.LowerBound[x.CurrentIterationCount]+LowerProblemTime
        x.History.LowerCount += 1
        #redirect_stdout(TT)
        PrintResults!(x,true)
      end
      x.History.CutCount[x.CurrentIterationCount] = x.CutIterations

      # checks for infeasibility stores solution
      if (x.CurrentLowerInfo.Feasibility)
        if (~x.ConvergenceCheck(x))

          # Solves upper bounding problem
          #redirect_stdout()
          #UpperProblemTime = @elapsed x.UpperProblem(x,CurrentNode)
          UpperProblemTime = x.UpperProblem(x,CurrentNode)
          #x.History.UpperTime[x.CurrentIterationCount] = x.History.UpperTime[x.CurrentIterationCount-1]+UpperProblemTime
          x.History.UpperCount += 1
          #redirect_stdout(TT)
          UBDTempSaveFlag = true
          PrintResults!(x,false)

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
          TT = stdout
          #redirect_stdout()
          PostprocessTime = x.Postprocess(x,CurrentNode) #@elapsed x.Postprocess(x,CurrentNode)
          #x.History.PostprocessTime[x.CurrentIterationCount] = x.History.PostprocessTime[x.CurrentIterationCount-1]+PostprocessTime
          PostSaveFlag = true
         # redirect_stdout(TT)

          # Checks to see if the node
          if (x.CurrentPostprocessInfo.Feasibility)
            if x.RepeatCheck(x,CurrentNode)
              SingleStorage!(x,CurrentNode)
            else
              Y1,Y2 = x.BisectionFunction(x,CurrentNode)
              x.NodeStorage(x,Y1,Y2,CurrentNode)
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
    (~UBDSaveFlag) && (x.History.UpperBound[x.CurrentIterationCount] = x.GlobalUpperBound[x.CurrentIterationCount-1])
    (~UBDTempSaveFlag) && (x.History.UpperTime[x.CurrentIterationCount] = x.History.UpperTime[x.CurrentIterationCount-1])
    (~PostSaveFlag) && (x.History.PostprocessTime[x.CurrentIterationCount] = x.History.PostprocessTime[x.CurrentIterationCount-1])
    x.History.Count[x.CurrentIterationCount] = x.CurrentNodeCount
    PrintIteration!(x)
    x.CurrentIterationCount += 1
  end

  x.SolutionValue = x.GlobalUpperBound           # Sets the solution value found
  PrintSolution!(x)                              # Prints the solution
end
