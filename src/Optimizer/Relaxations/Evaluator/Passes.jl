function SetValueConstruct(i::Int,N::Int,x_values::Vector{Float64},node::NodeBB)
    @inbounds MC{N}(x_values[i],x_values[i], Interval{Float64}(node.LowerVar[i],node.UpperVar[i]),seedg(Float64,i,N),seedg(Float64,i,N),false)
end

#=
function SetValuePost(x_values::Vector{Float64},val::MC{N},node::NodeBB) where N
    lower = val.cv
    upper = val.cc
    for i in 1:N
        @inbounds cv_val = cv_grad[i]
        @inbounds cc_val = cc_grad[i]
        @inbounds lower += cv_val > 0.0 ? cv_val*(node.LowerVar[i]-x_values[i]) : cv_val*(node.UpperVar[i]-x_values[i])
        @inbounds upper += cc_val > 0.0 ? cc_val*(node.UpperVar[i]-x_values[i]) : cc_val*(node.LowerVar[i]-x_values[i])
    end
    lower = max(lower,lo(val))
    upper = min(upper,hi(val))
    MC{N}(val.cv,val.cc,val.Interval(lower,upper),val.cv_grad,val.cc_grad,val.cnst)
end
=#
SetValuePost(x_values::Vector{Float64},val::MC{N},node::NodeBB) where N = val

function forward_eval(setstorage::Vector{T}, numberstorage::Vector{Float64}, setvalued::Vector{Bool},
                      nd::AbstractVector{JuMP.NodeData}, adj, const_values, parameter_values, current_node::NodeBB,
                      x_values::Vector{Float64}, subexpression_values, user_input_buffer = [];
                      user_operators::JuMP.Derivatives.UserOperatorRegistry=JuMP.Derivatives.UserOperatorRegistry()) where T

    @assert length(numberstorage) >= length(nd)
    @assert length(setstorage) >= length(nd)
    @assert length(setvalued) >= length(nd)

    children_arr = rowvals(adj)
    N = length(x_values)

    for k in length(nd):-1:1

        # compute the value of node k
        @inbounds nod = nd[k]
        if nod.nodetype == VARIABLE
            @inbounds setstorage[k] = SetValueConstruct(nod.index,N,x_values,current_node)
            setvalued[k] = true
        elseif nod.nodetype == VALUE
            @inbounds numberstorage[k] = const_values[nod.index]
            setvalued[k] = false
        elseif nod.nodetype == SUBEXPRESSION
            @inbounds isset = setvalued[nod.index]
            if isset
                @inbounds setstorage[k] = subexpression_values[nod.index]
            else
                @inbounds numberstorage[k] = subexpression_values[nod.index]
            end
            setvalued[k] == isset
        elseif nod.nodetype == PARAMETER
            @inbounds numberstorage[k] = parameter_values[nod.index]
            setvalued[k] = false
        elseif nod.nodetype == CALL
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            if op == 1 # :+
                tmp_sum = 0.0
                isset = true
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds tmp_sum += storage[ix]
                    @inbounds isset &= setvalued[ix]
                end
                setvalued[k] = isset
                if (isset)
                    setstorage[k] = tmp_sum
                else
                    numberstorage[k] = tmp_sum
                end
            elseif op == 2 # :-
                child1 = first(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                @inbounds isset = setvalued[ix1]
                @inbounds isset &= setvalued[ix2]
                @inbounds tmp_sub = storage[ix1]
                @inbounds tmp_sub -= storage[ix2]
                setvalued[k] = isset
                if (isset)
                    setstorage[k] = SetValuePost(x_values, tmp_sub, current_node)
                else
                    numberstorage[k] = tmp_sum
                end
            elseif op == 3 # :*
                tmp_prod = 1.0
                isset = true
                for c_idx in children_idx
                    @inbounds isset &= setvalued[children_arr[c_idx]]
                    if (isset)
                        @inbounds tmp_prod = SetValuePost(x_values, tmp_prod*storage[children_arr[c_idx]], current_node)
                    else
                        @inbounds tmp_prod *= storage[children_arr[c_idx]]
                    end
                end
                if (isset)
                    setstorage[k] = tmp_prod
                else
                    numberstorage[k] = tmp_prod
                end
            elseif op == 4 # :^
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                if (setvalued[ix1])
                    @inbounds base = setstorage[ix1]
                else
                    @inbounds exponent = numberstorage[ix1]
                end
                if (setvalued[ix2])
                    @inbounds base = setstorage[ix2]
                else
                    @inbounds exponent = numberstorage[ix2]
                end
                if exponent == 1
                    @inbounds storage[k] = base
                else
                    storage[k] = SetValuePost(x_values, pow(base,exponent), current_node)
                end
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds numerator = storage[ix1]
                @inbounds denominator = storage[ix2]
                storage[k] = SetValuePost(x_values, numerator/denominator, current_node)
            elseif op == 6 # ifelse
                @assert n_children == 3
                idx1 = first(children_idx)
                @inbounds condition = storage[children_arr[idx1]]
                @inbounds lhs = storage[children_arr[idx1+1]]
                @inbounds rhs = storage[children_arr[idx1+2]]
                storage[k] = SetValuePost(x_values, ifelse(condition == 1, lhs, rhs), current_node)
            elseif op >= USER_OPERATOR_ID_START
                evaluator = user_operators.multivariate_operator_evaluator[op - USER_OPERATOR_ID_START+1]
                f_input = view(user_input_buffer, 1:n_children)
                r = 1
                for c_idx in children_idx
                    ix = children_arr[c_idx]
                    f_input[r] = storage[ix]
                    r += 1
                end
                fval = MOI.eval_objective(evaluator, f_input)::T
                storage[k] = fval
                r = 1
                for c_idx in children_idx
                    ix = children_arr[c_idx]
                    r += 1
                end
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == CALLUNIVAR # univariate function
            op = nod.index
            @inbounds child_idx = children_arr[adj.colptr[k]]
            child_val = storage[child_idx]
            if op >= USER_UNIVAR_OPERATOR_ID_START
                userop = op - USER_UNIVAR_OPERATOR_ID_START + 1
                f = user_operators.univariate_operator_f[userop]
                fval = f(child_val)::T
            else
                fval = eval_univariate(op, child_val)
            end
            @inbounds storage[k] = SetValuePost(x_values, fval, current_node)
        elseif nod.nodetype == COMPARISON
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            result = true
            for r in 1:n_children-1
                cval_lhs = storage[children_arr[children_idx[r]]]
                cval_rhs = storage[children_arr[children_idx[r+1]]]
                if op == 1
                    result &= cval_lhs <= cval_rhs
                elseif op == 2
                    result &= cval_lhs == cval_rhs
                elseif op == 3
                    result &= cval_lhs >= cval_rhs
                elseif op == 4
                    result &= cval_lhs < cval_rhs
                elseif op == 5
                    result &= cval_lhs > cval_rhs
                end
            end
            storage[k] = SetValuePost(x_values, result, current_node)
        elseif nod.nodetype == LOGIC
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            cval_lhs = (storage[children_arr[first(children_idx)]] == 1)
            cval_rhs = (storage[children_arr[last(children_idx)]] == 1)
            if op == 1
                storage[k] = SetValuePost(x_values, cval_lhs && cval_rhs, current_node)
            elseif op == 2
                storage[k] = SetValuePost(x_values, cval_lhs || cval_rhs, current_node)
            end
        else
            error("Unrecognized node type $(nod.nodetype).")
        end
    end

    return storage[1]
end

function forward_eval_all(d::Evaluator,x)
    subexpr_values = d.subexpression_values
    user_operators = d.m.nlp_data.user_operators::JuMP.Derivatives.UserOperatorRegistry
    user_input_buffer = d.jac_storage
    for k in d.subexpression_order
        ex = d.subexpressions[k]
        subexpr_values[k] = forward_eval(ex.storage, ex.nd, ex.adj, ex.const_values,
                                         d.parameter_values, d.current_node,
                                         x, subexpr_values, user_input_buffer,
                                         user_operators=user_operators)
    end
    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.storage, ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x, subexpr_values, user_input_buffer,
                     user_operators=user_operators)
    end
    for ex in d.constraints
        forward_eval(ex.storage, ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x,subexpr_values, user_input_buffer,
                     user_operators=user_operators)
    end
end

function reverse_eval(reverse_storage::Vector{T},
                      nd::Vector{JuMP.NodeData}, adj) where T

    @assert length(reverse_storage) >= length(nd)
    @assert length(partials_storage) >= length(nd)

    reverse_storage[1] = one(T)

    for k in 2:length(nd)
        @inbounds nod = nd[k]
        if nod.nodetype == VALUE || nod.nodetype == LOGIC || nod.nodetype == COMPARISON || nod.nodetype == PARAMETER
            continue
        end
        @inbounds rev_parent = reverse_storage[nod.parent]
        @inbounds reverse_storage[k] = ifelse(rev_parent == 0.0, rev_parent, rev_parent*partial)
    end

    nothing
end

# looks good
function reverse_eval_all(d::Evaluator,x)
    for k in d.subexpression_order
        ex = d.subexpressions[k]
        reverse_eval(ex.reverse_storage, ex.nd, ex.adj)
    end
    if d.has_nlobj
        ex = d.objective
        reverse_eval(ex.reverse_storage, ex.nd, ex.adj)
    end
    for ex in d.constraints
        reverse_eval(ex.reverse_storage, ex.nd, ex.adj)
    end
    copyto!(d.last_x,x)
end

# looks good
function forward_reverse_pass(d::Evaluator,x)
    if d.last_x != x
        if d.has_reverse
            for i in d.fw_repeats
                forward_eval_all(d,x)
                reverse_eval_all(d,x)
            end
        else
            forward_eval_all(d,x)
        end
    end
end
