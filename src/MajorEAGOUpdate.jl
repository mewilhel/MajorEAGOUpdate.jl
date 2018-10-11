module MajorEAGOUpdate

import MathOptInterface

using Printf
using SparseArrays

using StaticArrays
using DiffRules
using Ipopt
using JuMP
using IntervalArithmetic

const MOI = MathOptInterface
const MOIU = MOI.Utilities

include("McCormick/McCormick.jl")
using .McCormick

include("Optimizer/IpoptSupplement.jl")
include("Optimizer/MathOptInterfaceEAGO.jl")

end # module
