using MathOptInterface
using MajorEAGOUpdate
using StaticArrays
using IntervalArithmetic
using Test
using JuMP
using Ipopt

const EAGO = MajorEAGOUpdate
const MOI = MathOptInterface

include("Optimizer/optimizer.jl")
include("McCormick/mccormick.jl")
