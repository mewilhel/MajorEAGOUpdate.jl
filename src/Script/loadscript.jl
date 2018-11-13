# Main decisions (replace nlp evaluator or recode the original EAGO/JuMP ones)

function loadscript!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function)
         # Adds a nonlinear constraint
         # Adds a nonlinear function
         # Initializes the NLP evaluator

end
