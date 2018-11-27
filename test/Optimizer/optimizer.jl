println("--- Test Branch and Bound Subroutines ---")
include("branch_bound.jl")

println("--- Test Domain Reduction Routines ---")
include("domain_reduction.jl")

println("--- Test MOI Interface Routines ---")
include("moi_interface.jl")

println("--- Test Other Utilities Routines ---")
include("other_utilities.jl")

println("--- Test Relaxation Routines ---")
include("relaxation_tools.jl")

#println("--- Test Main Solution Routines ---")
#include("optimize.jl")


#=
println("--- Test Implicit Optimization Routines ---")
include("implicit_tools.jl")
=#
