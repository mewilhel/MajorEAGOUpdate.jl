@testset "Test Base Implicit Routines" begin
    EAGO.MC_Differentiability!(0)
    opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)

    f(x,p) = x[1]*p[1]+p[1]
    g(x,p) = [x[1]*p[1]+p[1];
              x[1]*p[1]+2*p[1]]
    function h1(x,p)
        t1 = x[1]^2
        t2 = x[1]*p[1]
        t3 = 4.0
        t4 = t1 + t2
        t5 = t4 + t3
        return [t5]
    end

    hj1(x,p) = [2*x[1]+p[1]]

    P = [Interval(6.0,9.0)]
    X = [Interval(-0.78,-0.4)]
    p = [7.5]
    pmid = mid.(P)

    xIntv1 = Interval(1.0,3.0)
    xIBox = SVector{1,Interval{Float64}}([xIntv1])
    mBox = mid.(xIBox)
    np = 1
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    p_mc = [MC{np}(p[i],p[i],Interval(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]

    println("ran to 1")
    param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
    println("param: $param")

    println("ran to 2")
    hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
    println("hbnds: $hbnds")

    println("ran to 3")
    fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    println("fbnds: $fbnds")

    println("ran to 4")
    fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    println("fgbnds: $fgbnds")

    println("ran to 5")
    param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
    println("param: $param")

    println("ran to 6")
    @test isapprox(param[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(param[2][1].cv,-0.78,atol=1E-4)
    hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
    println("hbnds: $hbnds")

    @test isapprox(hbnds[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hbnds[1].cv,-0.6307841683146562,atol=1E-4)
    fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    println("fbnds: $fbnds")

    @test isapprox(fbnds.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fbnds.cv,2.5572882333553046,atol=1E-4)

    fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    println("fgbnds: $fbnds")
    @test isapprox(fgbnds[1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][2].cc,11.274523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][2].cv,10.057288233355305,atol=1E-4)
end
