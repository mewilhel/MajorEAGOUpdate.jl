module MajorEAGOUpdate

import MathOptInterface

using Printf
using SparseArrays

using JuMP
using JuMP.Derivatives

using Ipopt
using DiffRules
using Reexport
using StaticArrays
using IntervalArithmetic

const MOI = MathOptInterface
const MOIU = MOI.Utilities

include("McCormick/McCormick.jl")
using .McCormick

include("Optimizer/IpoptSupplement.jl")
include("Optimizer/MathOptInterfaceEAGO.jl")
include("Optimizer/Relaxations/Evaluator/Evaluator.jl")

end # module
