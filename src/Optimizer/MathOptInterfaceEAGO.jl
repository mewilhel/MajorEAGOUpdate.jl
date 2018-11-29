export NodeBB

# Access and manipulate the optimizer
include("MOIInterface/NodeBB.jl")
include("Relaxations/Scheme.jl")

include("MOIInterface/MOIInterface.jl")

# Static Storage and Problem Info Holders are set on optimize
include("Relaxations/RelaxModel.jl")

# Solution algorithm
include("BranchBound/BranchBound.jl")

# Domain Reduction Routines
include("DomainReduction/DBBT.jl")
include("DomainReduction/OBBT.jl")
include("DomainReduction/PoorLP.jl")
#include("DomainReduction/Quadratics.jl")

# Implicit Routines
include("Relaxations/ImplicitEvaluator/ImplicitEvaluator.jl")
