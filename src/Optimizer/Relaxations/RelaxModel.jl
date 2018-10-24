"""
    Takes an NLP data block structure, linear, quadratic constraints, and
    variable bounds and subsequently builds the relaxed model
"""
function RelaxModel!(src::Optimizer,trg,n::NodeBB,r::RelaxationScheme; load::Bool = false)
    if load
        # add linear terms to model
        RelaxLinear!(src,trg)
        # build NLP evaluator and save to EAGO object
        a = src.VariableNumber
        println("a: $a")
        b = src.NLPData
        println("b: $b")
        src.WorkingEvaluator = Build_NLP_Evaluator(MC{src.VariableNumber},src.NLPData.evaluator,src)
        # copy working evaluator into block if nonlinear block is needed
        if (r.OptimizerType == :NLP || r.OptimizerType == :MINLP)
            if ~isempty(src.NonlinearVariable)
                trg.nlp_data.evaluator = src.WorkingEvaluator
            end
        end
    else
        RelaxQuadratic!(src,trg,n,r,load)
        RelaxNonlinear!(src,trg,n,r,load)
        if ~isempty(src.NonlinearVariable)
            trg.nlp_data.evaluator.current_node = n
        end
    end
end
