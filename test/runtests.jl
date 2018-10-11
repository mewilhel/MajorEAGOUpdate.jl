using MathOptInterface
using MajorEAGOUpdate
using StaticArrays
using Test
using JuMP

const EAGO = MajorEAGOUpdate
const MOI = MathOptInterface

include("Optimizer/optimizer.jl")
include("McCormick/mccormick.jl")
