# OBBT
# Integer Variables: Fix one value and solve relaxation, fix other value and solve relaxtion, can dual simplex warm start potentially
# Apply to all nonlinear nonbinary variables
# Strong branching MINLP...

```
    Adds objective cut to constraint set
```
function SetObjectiveCut!(m::Optimizer)
    if typeof(m.Objective) == MOI.SingleVariable
        LoadObj = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm{Float64}(1.0,m.Objective.variable)],0.0)
    else
        LoadObj = m.Objective
    end
    m.ObjConstrIndx = MOI.ConstraintIndex[MOI.add_constraint(m.InitialRelaxedOptimizer, LoadObj, MOI.LessThan(m.GlobalUpperBound))]
end


```
    Excludes OBBT on variable indices that are tight for the solution of the relaxation (Should be done (for continuous)).
```
function TrivialFiltering!(x::Optimizer,y::NodeBB)
    TerminationStatus = MOI.get(x.WorkingRelaxedOptimizer, MOI.TerminationStatus())
    println("Trivial Term Status: $(TerminationStatus)")
    if (TerminationStatus != MOI.NoSolution)
        if (TerminationStatus == MOI.Success)
            StatPrimal = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
            println("Trivial Stat Primal: $(StatPrimal)")
            if StatPrimal == MOI.FeasiblePoint
                LowDelete = Int[]
                UppDelete = Int[]
                for i in x.OBBTWorkingLowerIndx
                    println("Lower Working Index: $(i)")
                    VarPrimal = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(),i)
                    println("Primal Variable: $(VarPrimal)")
                    (abs(VarPrimal[i] - y.LowerBnd[i]) <= x.OBBTTolerance) && push!(LowDelete,i)
                    println("Lower Delete: $(LowDelete)")
                end
                for i in x.OBBTWorkingUpperIndx
                    VarPrimal = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(),i)
                    println("Trivial Stat Primal: $(StatPrimal)")
                    (abs(y.UpperBnd[i] - VarPrimal[i]) <= x.OBBTTolerance) && push!(UppDelete,i)
                    println("UppDelete: $(UppDelete)")
                end
                deleteat!(x.OBBTWorkingLowerIndx,LowDelete)
                deleteat!(x.OBBTWorkingUpperIndx,UppDelete)
                println("End x.OBBTWorkingLowerIndx: $(x.OBBTWorkingLowerIndx)")
                println("End x.OBBTWorkingUpperIndx: $(x.OBBTWorkingUpperIndx)")
            end
        end
    end
end

#(Should be done (for continuous))
function AggressiveFiltering!(x::Optimizer)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    v = -ones(Float64,x.VariableNumber)

    # Copy prior index set (ignores linear and binary terms)
    OldLowIndx = x.OBBTWorkingLowerIndx
    OldUppIndx = x.OBBTWorkingUpperIndx
    NewLowIndx = x.OBBTWorkingLowerIndx
    NewUppIndx = x.OBBTWorkingUpperIndx

    # Exclude unbounded directions
    LowDelete = Int[]
    UppDelete = Int[]
    for i in NewLowIndx
        (y.LowerVar[i] == -Inf) && push!(LowDelete,i)
    end
    for i in NewUppIndx
        (y.UpperVar[i] == Inf) && push!(UppDelete,i)
    end
    deleteat!(NewLowIndx,LowDelete); LowDelete = Int[]
    deleteat!(NewUppIndx,UppDelete); UppDelete = Int[]

    # Begin the main algorithm
    for k in 1:x.OBBTAggrMaxIteration

        # Set index differences and vector for filtering direction
        LowIndxDiff = setdiff(OldLowIndx,NewLowIndx)
        UppIndxDiff = setdiff(OldUppIndx,NewUppIndx)
        LowDelete = Int[]
        UppDelete = Int[]
        for i in LowIndxDiff
            (v[i] < 0.0) && (v[i] = 0.0)
        end
        for i in UppIndxDiff
            (v[i] > 0.0) && (v[i] = 0.0)
        end

        # Termination Condition
        (isempty(union(NewLowIndx,NewUppIndx)) || (v == zeros(Float64,x.VariableNumber))) && break
        ((k >= 2) && (length(LowIndxDiff) + length(UppIndxDiff) < x.OBBTAggrMinDimLimit)) && break

        # Set objective in OBBT problem to filtering vector
        MOI.set(x.WorkingRelaxedOptimizer,
                 MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                 ScalarAffineFunction(ScalarAffineTerm.(v,x.OBBTvars),0.0))

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(x.WorkingRelaxedOptimizer)
        TerminationStatus = MOI.get(x.WorkingRelaxedOptimizer, MOI.TerminationStatus())
        if (TerminationStatus != MOI.NoSolution())
            if (TerminationStatus == MOI.Success())
                StatPrimal = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
                if StatPrimal == MOI.FeasiblePoint()
                    VarPrimal = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal())
                    NewLowIndx = copy(OldLowIndx)
                    NewUppIndx = copy(OldUppIndx)
                    for i in OldLowIndx
                        (VarPrimal[i] == y.LowerVar[i]) && push!(LowDelete,i)
                    end
                    for i in OldUppIndx
                        (VarPrimal[i] == y.UpperVar[i]) && push!(UppDelete,i)
                    end
                    deleteat!(NewLowIndx,LowDelete)
                    deleteat!(NewUppIndx,UppDelete)
                end
            else
                return false
            end
        else
            return false
        end
    end
    x.OBBTWorkingLowerIndx = NewLowIndx
    x.OBBTWorkingUpperIndx = NewUppIndx
    return true
end

```
    OBBT(x::EAGOOptimizer)

Performs OBBT with Filtering and Greedy Ordering
```
function OBBT(x::Optimizer,y::NodeBB)

    feasibility = true

    if ~x.FirstRelaxed
        RelaxModel!(x,y,x.WorkingRelaxedOptimizer)
    end

    # Sets indices to attempt OBBT on
    x.OBBTWorkingLowerIndx = copy(x.OBBTInitialLowerIndx)
    x.OBBTWorkingUpperIndx = copy(x.OBBTInitialUpperIndx)
    println("TestSet")
    println("x.OBBTInitialLowerIndx: $(x.OBBTInitialLowerIndx)")
    println("x.OBBTInitialUpperIndx: $(x.OBBTInitialUpperIndx)")
    println("x.OBBTWorkingLowerIndx: $(x.OBBTWorkingLowerIndx)")
    println("x.OBBTWorkingUpperIndx: $(x.OBBTWorkingUpperIndx)")

    # Set MINLP constraint that corresponds to objective
    for j in x.ObjConstrIndx
        println("Objective indx j: $(j)")
        MOI.set(x.WorkingRelaxedOptimizer, MOI.ConstraintSet(), j, MOI.LessThan(x.GlobalUpperBound))
    end

    # Prefiltering steps && and sets initial LP values
    TrivialFiltering!(x,y)
    println("Ran trivial filtering")
    if (x.OBBTAggrOn)
        feasibility = AggressiveFiltering!(x)
    end
    xLP = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(),x.VariableIndex)
    println("Got solution")

    while ~(isempty(x.OBBTWorkingLowerIndx) || isempty(x.OBBTWorkingUpperIndx))

        # Get lower value
        println("ran to me 1a")
        term1 = isempty(x.OBBTWorkingLowerIndx)
        println("ran to me 1b: xLP $xLP workingindex $(x.OBBTWorkingLowerIndx)")
        term2 = xLP[x.OBBTWorkingLowerIndx]
        println("ran to me 1c")
        term3 = y.LowerVar[x.OBBTWorkingLowerIndx]
        println("ran to me 1d")
        LowerValue = (isempty(x.OBBTWorkingLowerIndx)) ? Inf : minimum(xLP[x.OBBTWorkingLowerIndx] - y.LowerVar[x.OBBTWorkingLowerIndx])
        UpperValue = (isempty(x.OBBTWorkingUpperIndx)) ? Inf : minimum(y.UpperVar[x.OBBTWorkingUpperIndx] - xLP[x.OBBTWorkingUpperIndx])
        if (LowerValue <= UpperValue)

            # Get index
            i = minimum(xLP[x.OBBTWorkingLowerIndx] - LowerVar[x.OBBTWorkingLowerIndx])
            filter!(x == i, x.OBBTWorkingLowerIndx)

            # Solve optimization model
            MOI.set(x.WorkingRelaxedOptimizer, MOI.ObjectiveFunction{MOI.SingleVariable{Float64}}(), MOI.SingleVariable(x[i]))
            MOI.set(x.WorkingRelaxedOptimizer, ObjectiveSense(), MinSense)
            MOI.optimize!(x.WorkingOBBTOptimizer)
            TerminationStatus = MOI.get(x.WorkingRelaxedOptimizer, TerminationStatus())
            # If the solver succeeded check whether it was feasible, if infeasible stop algoritm
            if TerminationStatus == Success
                PrimalStatus = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
                if PrimalStatus == MOI.FeasiblePoint
                    xLP = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(),x.VariableIndex)
                    y.LowerVar[i] = IsIntegerVariable(x,i) ? ceil(xLP[i]) : xLP[i]
                elseif haskey(InfeasibleResultCodes,PrimalStatus)
                    feasibility = false
                    LowerBoundIndx = []
                    UpperBoundIndx = []
                end
            # If the solver results an infeasible code, check it to see whether it is guaranteed
            elseif haskey(InfeasibleTerminationCodes,TerminationStatus)
                if haskey(InfeasibleResultCodes,PrimalStatus)
                    feasibility = false
                    LowerBoundIndx = []
                    UpperBoundIndx = []
                end
            elseif haskey(ProblemTerminationCodes,TerminationStatus)
                if x.FlagSubSolverErrors
                    x.TerminationStatusCode = TerminationStatus
                    x.FailedSolver = true
                    error("The optimality-based bound tightening routine returned a
                           linear solve with the following status code: $(ErrorTerminationCodes[TerminationStatus]).
                           If this error is acceptable, set x.FlagSubSolverErrors = false to prevent an error from
                           being thrown.")
                end
            elseif haskey(ErrorTerminationCodes,TerminationStatus)
                x.TerminationStatusCode = TerminationStatus
                x.FailedSolver = true
                error("The optimality-based bound tightening routine returned a
                       linear solve with the following status code: $(ErrorTerminationCodes[TerminationStatus])")
            end
        else
            # Get index
            i = minimum(UpperVar[x.OBBTWorkingUpperIndx] - xLP[x.OBBTWorkingUpperIndx])
            filter!(x == i,  x.OBBTWorkingUpperIndx)

            # Solve optimization model
            set!(x.WorkingRelaxedOptimizer, ObjectiveSense(), MaxSense)
            MOI.set!(x.WorkingRelaxedOptimizer, MOI.ObjectiveFunction{MOI.SingleVariable{Float64}}(), MOI.SingleVariable(x[i]))
            MOI.optimize!(x.WorkingRelaxedOptimizer)
            TerminationStatus = MOI.get(x.WorkingRelaxedOptimizer, TerminationStatus())
            # If the solver succeeded check whether it was feasible, if infeasible stop algoritm
            if TerminationStatus == Success
                PrimalStatus = MOI.get(x.WorkingRelaxedOptimizer, MOI.PrimalStatus())
                if PrimalStatus == MOI.FeasiblePoint
                    xLP = MOI.get(x.WorkingRelaxedOptimizer, MOI.VariablePrimal(),x.VariableIndex)
                    UpperVar[i] = IsIntegerVariable(i) ? ceil(xLP[I]) : xLP[I]
                elseif haskey(InfeasibleResultCodes,PrimalStatus)
                    feasibility = false
                    LowerBoundIndx = []
                    UpperBoundIndx = []
                end
            # If the solver results an infeasible code, check it to see whether it is guaranteed
            elseif haskey(InfeasibleTerminationCodes,TerminationStatus)
                if haskey(InfeasibleResultCodes,PrimalStatus)
                    feasibility = false
                    LowerBoundIndx = []
                    UpperBoundIndx = []
                end
            elseif haskey(ProblemTerminationCodes,TerminationStatus)
                if x.FlagSubSolverErrors
                    x.TerminationStatusCode = TerminationStatus
                    x.FailedSolver = true
                    error("The optimality-based bound tightening routine returned a
                           linear solve with the following status code: $(ErrorTerminationCodes[TerminationStatus]).
                           If this error is acceptable, set x.FlagSubSolverErrors = false to prevent an error from
                           being thrown.")
                end
            elseif haskey(ErrorTerminationCodes,TerminationStatus)
                x.TerminationStatusCode = TerminationStatus
                x.FailedSolver = true
                error("The optimality-based bound tightening routine returned a
                       linear solve with the following status code: $(ErrorTerminationCodes[TerminationStatus])")
            end
        end
        #GenerateLVB!(x)
        TrivialFiltering!(x,y)
        end
    return feasibility
end

```

Add an objective function cutting constraint to the model
```
function AddObjectiveCut!(m::Optimizer)
    if (typeof(m.Objective) == MOI.SingleVariable)
        VarToAffine = MOI.ScalarAffineTerm{Float64}(1.0,m.Objective.variable)
        m.ObjConstrIndx = MOI.ConstraintIndex[MOI.add_constraint(m.InitialRelaxedOptimizer, VarToAffine, MOI.LessThan(m.GlobalUpperBound))]
    elseif ((typeof(m.Objective) == MOI.ScalarQuadraticFunction{Float64}) || (typeof(m.Objective) == MOI.ScalarAffineFunction{Float64}))
        m.ObjConstrIndx = MOI.ConstraintIndex[MOI.add_constraint(m.InitialRelaxedOptimizer, m.Objective, MOI.LessThan(m.GlobalUpperBound))]
    end
end



    #=
    LimitTerminationCodes = Dict{MOI.TerminationStatusCode,Symbol}()          # The solver hit a limit
    InfeasibleTerminationCodes = Dict{MOI.TerminationStatusCode,Symbol}()     # The result is infeasible
    ProblemTerminationCodes = Dict{MOI.TerminationStatusCode,Symbol}()        # There was a numerical issue in subproblem
    ErrorTerminationCodes = Dict{MOI.TerminationStatusCode,Symbol}()          # An error occured that should terminate the algorithm

    LimitTerminationCodes[MOI.IterationLimit] = :IterationLimit
    LimitTerminationCodes[MOI.TimeLimit] = :TimeLimit
    LimitTerminationCodes[MOI.NodeLimit] = :NodeLimit
    LimitTerminationCodes[MOI.SolutionLimit] = :SolutionLimit
    LimitTerminationCodes[MOI.MemoryLimit] = :MemoryLimit
    LimitTerminationCodes[MOI.NormLimit] = :NormLimit
    LimitTerminationCodes[MOI.OtherLimit] = :OtherLimit

    InfeasibleTerminationCodes[MOI.InfeasibleNoResult] = :InfeasibleNoResult
    InfeasibleTerminationCodes[MOI.InfeasibleOrUnbounded] = :InfeasibleOrUnbounded

    ProblemTerminationCodes[MOI.UnboundedNoResult] = :UnboundedNoResult
    ProblemTerminationCodes[MOI.NumericalError] = :NumericalError
    ProblemTerminationCodes[MOI.AlmostSuccess] = :AlmostSuccess
    ProblemTerminationCodes[MOI.SlowProgress] = :SlowProgress

    ErrorTerminationCodes[MOI.InvalidModel] = :InvalidModel
    ErrorTerminationCodes[MOI.InvalidOption] = :InvalidOption
    ErrorTerminationCodes[MOI.OtherError] = :OtherError
    ErrorTerminationCodes[MOI.Interrupted] = :Interrupted

    InfeasibleResultCodes = Dict{MOI.ResultStatusCode,Symbol}()
    UnknownResultCodes = Dict{MOI.ResultStatusCode,Symbol}()

    InfeasibleResultCodes[MOI.InfeasiblePoint] = :InfeasiblePoint
    InfeasibleResultCodes[MOI.InfeasibilityCertificate] = :InfeasibilityCertificate

    UnknownResultCodes[MOI.NearlyFeasiblePoint] = :NearlyFeasiblePoint
    UnknownResultCodes[MOI.NearlyInfeasibilityCertificate] = :NearlyInfeasibilityCertificate
    UnknownResultCodes[MOI.ReductionCertificate] = :ReductionCertificate
    UnknownResultCodes[MOI.NearlyReductionCertificate] = :NearlyReductionCertificate
    UnknownResultCodes[MOI.UnknownResultStatus] = :UnknownResultStatus
    UnknownResultCodes[MOI.OtherResultStatus] = :OtherResultStatus
    =#
