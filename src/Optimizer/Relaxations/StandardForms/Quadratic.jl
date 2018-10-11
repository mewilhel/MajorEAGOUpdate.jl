function RelaxQuadratic!(src::Optimizer,trg::T,n::NodeBB) where {T<:MOI.AbstractOptimizer}
    midn = (n.UpperVar - n.LowerVar)/2.0
    # Relax Convex Quadratic Terms
    for (func,set,ind) in src.QuadraticLEQConstraints

        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.QuadraticGEQConstraints
        MOI.add_constraint(trg, func, set)
    end
    for (func,set,ind) in src.QuadraticEQConstraints
        MOI.add_constraint(trg, func, set)
    end
    # Relax Remainging Terms via McCormick Inequalities
end

mutable struct ScalarQuadraticFunction{T} <: AbstractScalarFunction
    affine_terms::Vector{ScalarAffineTerm{T}}
    quadratic_terms::Vector{ScalarQuadraticTerm{T}}
    constant::T
end

struct ScalarQuadraticTerm{T}
    coefficient::T
    variable_index_1::VariableIndex
    variable_index_2::VariableIndex
end

struct ScalarAffineTerm{T}
    coefficient::T
    variable_index::VariableIndex
end

temp = ScalarAffineFunction{Float64}()
ScalarAffineFunction{Float64}.terms = Vector{ScalarAffineTerm{T}}

constant = func.constant
for qt in quadratic_terms
    constant += qt.coefficient*midn[VItoSto[qt.variable_index_1]]*midn[VItoSto[qt.variable_index_2]]
end
ScalarAffineFunction{Float64}(affine_terms, constant)
