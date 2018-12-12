"""
    Takes an NLP data block structure, linear, quadratic constraints, and
    variable bounds and subsequently builds the relaxed model
"""
function RelaxModel!(src::Optimizer,trg,n::NodeBB,r::RelaxationScheme; load::Bool = false)
    if load
        # add linear terms to model
        RelaxLinear!(src,trg)
        if (typeof(src.NLPData.evaluator) != MajorEAGOUpdate.EmptyNLPEvaluator)
            # copy working evaluator into block if nonlinear block is needed
            if (r.OptimizerType == :NLP || r.OptimizerType == :MINLP)
                if ~isempty(src.NonlinearVariable)
                    trg.nlp_data = src.WorkingEvaluatorBlock
                end
            end
        end
    else
        ~isinf(src.GlobalUpperBound) && ObjCutLinear!(src,trg)
        RelaxQuadratic!(trg,src,n,r)
        if ~isempty(src.NonlinearVariable)
            if MOI.supports(trg, MOI.NLPBlock())
                evaluator = src.WorkingEvaluatorBlock.evaluator
                set_current_node!(evaluator,n)
                nlp_data = MOI.NLPBlockData(src.NLPData.constraint_bounds,
                                            evaluator,
                                            src.NLPData.has_objective)
                MOI.set(trg, MOI.NLPBlock(), nlp_data)
            else
                MidPointAffine!(src,trg,n,r)
            end
        end
    end
end
