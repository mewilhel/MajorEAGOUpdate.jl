X = EAGO.Optimizer()
sip1 = SIP_opts(X)

# solves example SIP #1 with DAG contractor disabled
    EAGO = MajorEAGOUpdate
    SIPopt1 = SIP_opts()

    function f1(x)
        (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    end
    function gSIP1(x,p)
        (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
    end
    Xlow = [-1000.0,-1000.]
    Xhigh = [1000.0,1000.0]
    Plow = [0.0]
    Phigh = [1.0]
    SIPoutput1 = Explicit_SIP_Solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, SIPopt1)

    SIPopt1 = SIP_opts()

    f1(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    gSIP1(x,p) = (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
    Xlow = [-1000.0,-1000.0]
    Xhigh = [1000.0,1000.0]
    Plow = [0.0]
    Phigh = [1.0]

    println("error test 1")
    SIPopt1.r0 = 2.0
    SIPopt1.eps_g0 = -0.9
    EAGO.Explicit_SIP_Solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, SIPopt1)

    println("error test 2")
    SIPopt1.r0 = 0.1
    SIPopt1.eps_g0 = 0.9
    EAGO.Explicit_SIP_Solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, SIPopt1)

    println("error test 3")
    SIPopt1.r0 = 2.0
    SIPopt1.LBP_Opt.DAG_depth = 100
    EAGO.Explicit_SIP_Solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, SIPopt1)
