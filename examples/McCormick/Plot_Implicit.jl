using MajorEAGOUpdate

const EAGO = MajorEAGOUpdate

using DataFrames
using IntervalArithmetic
using StaticArrays
using CSV
using Plots
#using PyPlot
#pyplot()

opts1 =  mc_opts(0.5,2,:Dense,:Newton,1,1,0.0)
function h1(x,p)
    return [x[1]^2+x[1]*p[1]+4.0]
end
hj1(x,p) = [2.0*x[1]+p[1]]

pl = 6.0; pu = 9.0; xl = -0.78; xu = -0.40
P = [Interval(pl,pu)]
pmid = [7.5]
intv = 0.01
np = Integer((pu-pl)/intv+1)
X = [Interval(xl,xu)]

x_grid = zeros(Float64,np)
vvalues = zeros(Float64,np)
solbranch = zeros(Float64,np)
solbranch2 = zeros(Float64,np)

sto_cc_grad = [0.0]
sto_cv_grad = [0.0]

sto_cc_grada = [0.0]
sto_cv_grada = [0.0]

ccvalNSnr = zeros(Float64,np)  # Storage for nonsmooth, subgradient refinement
cvvalNSnr = zeros(Float64,np)
ccvalNSnn = zeros(Float64,np)  # Storage for nonsmooth, no subgradient refinement
cvvalNSnn = zeros(Float64,np)

ccvalNSnr2 = zeros(Float64,np)  # Storage for nonsmooth, subgradient refinement
cvvalNSnr2= zeros(Float64,np)
ccvalNSnn2 = zeros(Float64,np)  # Storage for nonsmooth, no subgradient refinement
cvvalNSnn2 = zeros(Float64,np)

tempbox = SVector{1,Interval{Float64}}(1.0*[Interval(xl,xu)])
tempmid =  mid.(tempbox)

EAGO.MC_Differentiability!(0)
param_ns_sr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
temp1 = [MC{1}(pmid[1],pmid[1],Interval(pl,pu),seedg(Float64,1,1),seedg(Float64,1,1), false)]
temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
#=
for (count1,j) in enumerate(pl:intv:pu)
    println("count1,j: $count1, $j")
    vvalues[count1] = j
    solbranch[count1] = (j-sqrt(j^2-16.0))/2.0
    solbranch2[count1] = (-j+sqrt(j^2-16.0))/2.0
    temp1 = [MC{1}(Float64(j),Float64(j),Interval(pl,pu),seedg(Float64,1,1),seedg(Float64,1,1),false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
    ccvalNSnn[count1] = temp1_MC[end].cc
    cvvalNSnn[count1] = temp1_MC[end].cv
end
=#
#=
param_ns_nr1 = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
param_ns_nr2 = GenExpansionParams(h2,hj2,X,P,pmid,opts1)
for (count1,j) in enumerate(xl:intv:xu)
    temp1 = [MC{1}(Float64(j),Float64(j),Interval(pl,pu),seedg(Float64,1,1),seedg(Float64,1,1),false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_nr1)
    temp2_MC = MC_impRelax(h2,hj2,temp1,pmid,X,P,opts1,param_ns_nr2)
    if (count1 == 31)
        println("temp1 at 31: $(temp1_MC)")
        sto_cc_grad[1] = temp1_MC[1].cc_grad[1]
        sto_cv_grad[1] = temp1_MC[1].cv_grad[1]
        sto_cc_grada[1] = temp2_MC[1].cc_grad[1]
        sto_cv_grada[1] = temp2_MC[1].cv_grad[1]
    end
    ccvalNSnr[count1] = temp1_MC[end].cc
    cvvalNSnr[count1] = temp1_MC[end].cv
    ccvalNSnr2[count1] = temp2_MC[end].cc
    cvvalNSnr2[count1] = temp2_MC[end].cv
end
=#
