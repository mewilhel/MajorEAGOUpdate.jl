function EAGODefault_PreProcess!(x::Optimizer,y::NodeBB)
    # Sets initial feasibility
    feas = true; rept = 0

    # Turn off tightening if too deep
    LPFlag = (x.PoorManLPDepth > x.CurrentIterationCount) ? false : true
    UQFlag = (x.UniQuadDepth > x.CurrentIterationCount) ? false : true
    BQFlag = (x.BiQuadDepth > x.CurrentIterationCount) ? false : true
    CPFlag = (x.CPWalkDepth > x.CurrentIterationCount) ? false : true
    x.OBBTActiveCurrent = (x.OBBTDepth > x.CurrentIterationCount) ? false : true
    ContFlag = (LPFlag || UQFlag || BQFlag || CPFlag)
    #=
    while ContFlag
        # Turns off tightening if number of reptitions exceeded
        (x.PoorManLPRepts > rept) && (LPFlag = false)
        (x.UniQuadRepts > rept) && (UQFlag = false)
        (x.BiQuadRepts > rept) && (BQFlag = false)
        (x.CPWalkRepts > rept) && (CPFlag = false)

        # Runs tightening routines, stops if infeasibility proven
        #LPFlag && (feas = PoorManLP(x,y));                (~feas) && (break)
        #UQFlag && (feas = UnivariateQuadratic(x,y));      (~feas) && (break)
        #BQFlag && (feas = BivariateQuadratic(x,y));       (~feas) && (break)
        #CPFlag && (feas = CPWalk!(x,y));               (~feas) && (break)

        ContFlag = (LPFlag || UQFlag || BQFlag || CPFlag)
        rept += 1
    end
    =#

    #=
    rept = 0
    if x.OBBTActiveCurrent
        if x.InitialOBBTOptimizer != DummyOptimizer()
            MOI.copy!(x.InitialOBBTOptimizer,x.WorkingOBBTOptimizer)
        end
        RelaxModel!(x,y,x.WorkingOBBTOptimizer)
    end
    =#
    #=
    while x.OBBTActiveCurrent
        (x.OBBTRepts > rept) && (x.OBBTActiveCurrent = false)
        x.OBBTActiveCurrent && (feas = OBBT(x,y));            (~feas) && (break)
        rept += 1
    end
=#
    x.CurrentPreprocessInfo.Feasibility = feas
end
