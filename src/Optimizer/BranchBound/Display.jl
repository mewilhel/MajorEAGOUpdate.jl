"""
    EAGO.print_sol!(x::EAGOOptimizer)

Prints solution information for the B&B problem. Displays first node found, solution value,
solution, and time spent solving subproblems.
"""
function PrintSolution!(x::Optimizer)
  temp::Float64 = 0.0
  println("First Solution Found at Node $(x.FirstSolutionNode)")
  if (x.Verbosity == 1 || x.Verbosity == 2 || x.Verbosity == 3)
    println("UBD = $(x.GlobalUpperBound)")
    println("Solution is :")
    if (x.FeasibleSolutionFnd)
      for i=1:length(x.ContinuousSolution)
        temp = x.ContinuousSolution[i]
        println("    X[$i] = $temp")
      end
    end
    #println("Total LBD problems solved = $(x.History.LowerCount) in $(x.History.LowerTime[x.CurrentIterationCount]) seconds.")
  #  println("Total UBD problems solved = $(x.History.UpperCount) in $(x.History.UpperTime[x.CurrentIterationCount]) seconds.")
  #  println("Total time spent preprocessing =  $(x.History.PreprocessTime[x.CurrentIterationCount]) seconds.")
  #  println("Total time spent postprocessing = $(x.History.PostprocessTime[x.CurrentIterationCount]) seconds.")
  end
end

"""
    EAGO.print_node!(x::BnBSolver,id::Int64,lbd::Float64,box::Vector{Interval{V}}) where {V<:AbstractFloat}

Prints node information for the B&B problem. Node id, bound, and interval box.
"""
function PrintNode!(id::Int,x::NodeBB)
    println("Node ID: $(id), Lower Bound: $(x.LowerBound), Lower Variable Bounds: $(x.LowerVar), Upper Variable Bounds: $(x.UpperVar)")
end


"""
    EAGO.print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64,ubd::Float64,feasL::Bool,feasU::Bool)

Prints the iteration information if the Verbosity is set to "Normal" or "Full".
The header is displayed every hdr_intv, the iteration info is displayed every
itr_intv
"""
function PrintIteration!(x::Optimizer)
  if (x.Verbosity == 1 || x.Verbosity == 2 || x.Verbosity == 3)
    # prints header line every B.hdr_intv times
    if (mod(x.CurrentIterationCount,x.HeaderInterations) == 0 || x.CurrentIterationCount == 1)
      println("Iteration   NodeID    Current_LBD     Global_LBD     Global_UBD      NodesLeft     Absolute_Gap    Absolute_Ratio     LBD_Feas     UBD_Feas")
    end
    # prints iteration summary every B.itr_intv times
    sbool1 = x.CurrentLowerInfo.Feasibility ? "true" : "false"
    sbool2 = x.CurrentUpperInfo.Feasibility ? "true" : "false"
    nid = x.GlobalUpperBound
    lbdp = x.CurrentLowerInfo.Feasibility ? x.CurrentLowerInfo.Value : Inf
    if ((mod(x.CurrentIterationCount,x.OutputInterations) == 0))
      ptr_arr1 = join([Printf.@sprintf("%6u",x) for x in Int[x.CurrentIterationCount x.CurrentNodeCount]], ",   ")
      ptr_arr2 = join([Printf.@sprintf("%3.7f",x) for x in Float64[lbdp x.GlobalLowerBound x.GlobalUpperBound]], ",     ")
      ptr_arr3 = join([Printf.@sprintf("%6u",x) for x in Int[x.CurrentNodeCount]], ",")
      ptr_arr4 = join([Printf.@sprintf("%3.7f",x) for x in Float64[abs(x.GlobalUpperBound-x.GlobalLowerBound) abs(x.GlobalUpperBound-x.GlobalLowerBound)/(min(abs(x.GlobalLowerBound),abs(x.GlobalUpperBound)))]], ",       ")
      ptr_arr5 = join([Printf.@sprintf("%s",x) for x in String[sbool1 sbool2]], ",       ")
#      ptr_arr1 = join([@sprintf("%6u",x) for x in Int[k_int nid]], ",   ")
#      ptr_arr2 = join([@sprintf("%3.7f",x) for x in Float64[lbdp lbd ubd]], ",     ")
#      ptr_arr3 = join([@sprintf("%6u",x) for x in Int[k_nod]], ",")
#      ptr_arr4 = join([@sprintf("%3.7f",x) for x in Float64[abs(ubd-lbd) abs(ubd-lbd)/(min(abs(lbd),abs(ubd)))]], ",       ")
#      ptr_arr5 = join([@sprintf("%s",x) for x in Bool[sbool1 sbool2]], ",       ")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3),",        ", ptr_arr4,",        ", ptr_arr5)
    end
  end
end

"""
    EAGO.PrintResults!(B::EAGOOptimizer,lbd_bool::Bool)

Prints the results of a single bounding problem.
"""
function PrintResults!(B::Optimizer,lbd_bool::Bool)
  if (B.Verbosity == 2 || B.Verbosity == 3)
    if (lbd_bool)
      println("Lower Bound: $(B.CurrentLowerInfo.Value), Solution: $(B.CurrentLowerInfo.Solution), Feasibility: $(B.CurrentLowerInfo.Feasibility)")
    else
      println("Upper Bound: $(B.CurrentUpperInfo.Value), Solution: $(B.CurrentUpperInfo.Solution), Feasibility: $(B.CurrentUpperInfo.Feasibility)")
    end
  end
end
