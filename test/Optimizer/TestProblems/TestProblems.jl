# Run LP tests
#=
# One to four pass
for i in 1:4
    include("LP/Prob$i.jl")
end
=#

# Run QP tests
#include("QP/Prob1.jl")
#include("QP/Prob3.jl")
#=
# Problem 1 ---- 2 Passes on NL
# Problem 3 Fails on both NL and quad handling
for i in 1:2
    include("QP/Prob$i.jl")
end
=#

# Run NLP tests
include("NLP/Prob7.jl")
#=
# Problem 1-3 Passes
for i in 1:3
    include("NLP/Prob$i.jl")
end
=#
# Run Implicit test problems
#include("Implicit/Ex5_1.jl")
#include("Ex5_2.jl")
