"""
    Takes an NLP data block structure, linear, quadratic constraints, and
    variable bounds and subsequently builds the relaxed model
"""
function RelaxModel!(src::Optimizer,trg::T,n::NodeBB,r::RelaxationScheme; load::Bool = false) where {T<:MOI.AbstractOptimizer}
    if load
        RelaxLinear!(src,trg)NLPData
        trg.d = Build_NLP_Evaluator(src.NLPData.evaluator)
    else
        RelaxQuadratic!(src,trg,n,r,load)
        RelaxNonlinear!(src,trg,n,r,load)
    end
end



function RelaxNonlinear!(src::Optimizer,trg::T,r::RelaxationScheme, load::Bool) where {T<:MOI.AbstractOptimizer}
    for i in NonlinearRelaxedAfterLoad
    end
end

function RelaxLinear(m::Optimizer)
end
# DefaultRelaxationScheme()
