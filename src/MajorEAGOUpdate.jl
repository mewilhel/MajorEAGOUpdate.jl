module MajorEAGOUpdate

import MathOptInterface

using Printf
using SparseArrays
using LinearAlgebra

using JuMP
using JuMP.Derivatives

using Ipopt
using Clp
using DiffRules
using Reexport
using StaticArrays
using IntervalArithmetic
using Calculus

const MOI = MathOptInterface
const MOIU = MOI.Utilities

include("McCormick/McCormick.jl")
using .McCormick

import Base.eltype

include("Optimizer/IpoptSupplement.jl")
include("Optimizer/MathOptInterfaceEAGO.jl")

include("Optimizer/Relaxations/StandardForms/Linear.jl")
include("Optimizer/Relaxations/StandardForms/Quadratic.jl")
include("Optimizer/Relaxations/StandardEvaluator/Evaluator.jl")

export SIP_opts, Explicit_SIP_Solve

include("SemiInfinite/SemiInfinite.jl")

end # module
