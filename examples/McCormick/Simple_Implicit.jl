using MajorEAGOUpdate
using IntervalArithmetic
using StaticArrays

const EAGO = MajorEAGOUpdate

EAGO.MC_Differentiability!(0)
opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)
#=
generates the expansion point parameters for the function using the opts
options using inverse preconditioner
=#
f(x,p) = x[1]*p[1]+p[1]
g(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h1(x,p)
    #println("(h1) x: $x, p: $p")
    #println("p: $p")
    t1 = x[1]^2
    #println("t1 $t1")
    t2 = x[1]*p[1]
    #println("t2: $t2")
    t3 = 4.0
    t4 = t1 + t2
    #println("t4: $t4")
    t5 = t4 + t3
    #println("t5: $t5")
    return [t5]
end
function hj1(x,p)
    #println("(hj1) x: $x, p: $p")
    t6 = 2*x[1]+p[1]
    #println("t6: $t6")
    [t6]
end
P = [Interval(6.0,9.0)]
X = [Interval(-0.78,-0.4)]
p = [7.5]
pmid = mid.(P)

#=
relaxes the equality h(x,p)
=#
np = 1
szero = @SVector zeros(np)
sone = @SVector ones(np)
nxi = 1
p_mc = [MC{np}(p[i],p[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]
xa_mc = [(MC{np}(X[i].lo,X[i].lo,Interval{Float64}(X[i].lo,X[i].lo),szero,szero,false)) for i=1:nxi]
xA_mc = [(MC{np}(X[i].hi,X[i].hi,Interval{Float64}(X[i].hi,X[i].hi),szero,szero,false)) for i=1:nxi]
z_mc = 0.5*xa_mc+0.5*xA_mc

param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
