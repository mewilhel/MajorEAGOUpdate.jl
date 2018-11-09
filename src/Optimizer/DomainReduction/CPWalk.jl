# Performs constraint walking on nonlinear terms
function CPWalk!(x::Optimizer,n::NodeBB)

    # runs at midpoint bound
    midx = n.LowerVar + (n.UpperVar - n.LowerVar)/2.0

    # set working node to n, copies pass parameters from EAGO optimizer
    x.WorkingEvaluatorBlock.evaluator.current_node = n
    

    # Run forward-reverse pass and retreive node
    forward_reverse_pass(x.WorkingEvaluatorBlock.evaluator, midx)
    n = x.WorkingEvaluatorBlock.current_node

    # if interval bounds empty then label as infeasible

end
