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

function forward_eval(setstorage::Vector{T}, numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                      nd::AbstractVector{JuMP.NodeData}, adj, const_values, parameter_values, current_node::NodeBB,
                      x_values::Vector{Float64}, subexpression_values, user_input_buffer;
                      user_operators::JuMP.Derivatives.UserOperatorRegistry=JuMP.Derivatives.UserOperatorRegistry()) where T

    @assert length(numberstorage) >= length(nd)
    @assert length(setstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    children_arr = rowvals(adj)
    N = length(x_values)

    for k in length(nd):-1:1
        # compute the value of node k
        @inbounds nod = nd[k]
        if nod.nodetype == JuMP.Derivatives.VARIABLE
            setstorage[k] = SetValueConstruct(nod.index,N,x_values,current_node)
            println("AT k: $k,   VARIABLE assigment yeilds setstorage: $(setstorage[k])")
            numvalued[k] = false
        elseif nod.nodetype == JuMP.Derivatives.VALUE
            @inbounds numberstorage[k] = const_values[nod.index]
            println("AT k: $k,   Value assignment yeilds numstorage: $(numberstorage[k])")
            numvalued[k] = true
        elseif nod.nodetype == JuMP.Derivatives.SUBEXPRESSION
            @inbounds isnum = numvalued[nod.index]
            if isnum
                @inbounds numberstorage[k] = subexpression_values[nod.index]
            else
                @inbounds setstorage[k] = subexpression_values[nod.index]
            end
            println("AT k: $k,   Subexpression evaluation yeilds subexpr: $(subexpression_values[nod.index])")
            numvalued[k] == isnum
        elseif nod.nodetype == JuMP.Derivatives.PARAMETER
            @inbounds numberstorage[k] = parameter_values[nod.index]
            println("AT k: $k,   Parameter assignment yeilds numstorage: $(numberstorage[k])")
            numvalued[k] = true
        elseif nod.nodetype == JuMP.Derivatives.CALL
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            if op == 1 # :+
                tmp_sum = 0.0
                isnum = true
                chdset = true
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds chdset = numvalued[c_ix]
                    if (chdset)
                        @inbounds tmp_sum += setstorage[c_ix]
                    else
                        @inbounds tmp_sum += numberstorage[c_ix]
                    end
                    @inbounds isnum &= chdset
                end
                numvalued[k] = isnum
                if (isnum)
                    numberstorage[k] = tmp_sum
                else
                    setstorage[k] = tmp_sum
                end
                println("AT k: $k,   PLUS CALL assigment yeilds setstorage: $(tmp_sum)")
            elseif op == 2 # :-
                child1 = first(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                @inbounds isnum = chdset1
                @inbounds isnum &= chdset2
                if chdset1
                    @inbounds tmp_sub = numberstorage[ix1]
                else
                    @inbounds tmp_sub = setstorage[ix1]
                end
                if chdset2
                    @inbounds tmp_sub -= numberstorage[ix2]
                else
                    @inbounds tmp_sub -= setstorage[ix2]
                end
                numvalued[k] = isnum
                if (isnum)
                    numberstorage[k] = tmp_sum
                    println("AT k: $k,   MINUS CALL assigment yeilds setstorage: $(tmp_sum)")
                else
                    setstorage[k] = SetValuePost(x_values, tmp_sub, current_node)
                    println("AT k: $k,   MINUS CALL assigment yeilds setstorage: $(setstorage[k])")
                end
            elseif op == 3 # :*
                tmp_prod = 1.0
                isnum = true
                chdset = true
                for c_idx in children_idx
                    chdset = numvalued[children_arr[c_idx]]
                    println("references: $(children_arr[c_idx])")
                    println("chdset: $(chdset)")
                    isnum &= chdset
                    if (chdset)
                        tmp_prod *= numberstorage[children_arr[c_idx]]
                        println("tmp_prod: $tmp_prod")
                    else
                        tmp_prod = tmp_prod*setstorage[children_arr[c_idx]]
                        println("tmp_prod: $tmp_prod")
                    end
                end
                if (isnum)
                    numberstorage[k] = tmp_prod
                    println("AT k: $k,   MULT CALL assigment yeilds numberstorage: $(tmp_prod)")
                else
                    setstorage[k] = SetValuePost(x_values, tmp_prod, current_node)
                    println("AT k: $k,   MULT CALL assigment yeilds setstorage: $(tmp_prod)")
                end
                numvalued[k] = isnum
            elseif op == 4 # :^
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    @inbounds base = numberstorage[ix1]
                else
                    @inbounds base = setstorage[ix1]
                end
                if chdset2
                    @inbounds exponent = numberstorage[ix2]
                else
                    @inbounds exponent = setstorage[ix2]
                end
                if exponent == 1
                    if chdset1
                        @inbounds numberstorage[k] = base
                    else
                        @inbounds setstorage[k] = base
                    end
                else
                    if chdset1 || chdset2
                        setstorage[k] = SetValuePost(x_values, pow(base,exponent), current_node)
                        println("AT k: $k,   EXP CALL assigment yeilds setstorage: $(setstorage[k])")
                    else
                        numberstorage[k] = pow(base,exponent)
                        println("AT k: $k,   EXP CALL assigment yeilds setstorage: $(numberstorage[k])")
                    end
                end
                numvalued[k] = ~(chdset1 || chdset2)
            elseif op == 5 # :/
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    @inbounds numerator = setstorage[ix1]
                else
                    @inbounds numerator = numberstorage[ix1]
                end
                if chdset2
                    @inbounds denominator = setstorage[ix2]
                else
                    @inbounds denominator = numberstorage[ix2]
                end
                if chdset1 || chdset2
                    numberstorage[k] = numerator/denominator
                    println("AT k: $k,   DIV CALL assigment yeilds setstorage: $(numberstorage[k])")
                else
                    setstorage[k] = SetValuePost(x_values, numerator/denominator, current_node)
                    println("AT k: $k,   DIV CALL assigment yeilds setstorage: $(setstorage[k])")
                end
                numvalued[k] = chdset1 || chdset2
            elseif op == 6 # ifelse
                @assert n_children == 3
                idx1 = first(children_idx)
                @inbounds chdset1 = numvalued[idx1]
                if chdset1
                    @inbounds condition = setstorage[children_arr[idx1]]
                else
                    @inbounds condition = numberstorage[children_arr[idx1]]
                end
                @inbounds chdset2 = numvalued[children_arr[idx1+1]]
                @inbounds chdset3 = numvalued[children_arr[idx1+2]]
                if chdset2
                    @inbounds lhs = setstorage[children_arr[idx1+1]]
                else
                    @inbounds lhs = numberstorage[children_arr[idx1+1]]
                end
                if chdset3
                    @inbounds rhs = setstorage[children_arr[idx1+2]]
                else
                    @inbounds rhs = numberstorage[children_arr[idx1+2]]
                end
                error("IF ELSE TO DO")
                #storage[k] = SetValuePost(x_values, ifelse(condition == 1, lhs, rhs), current_node)
            elseif op >= JuMP.Derivatives.USER_OPERATOR_ID_START
                evaluator = user_operators.multivariate_operator_evaluator[op - USER_OPERATOR_ID_START+1]
                f_input = view(user_input_buffer, 1:n_children)
                r = 1
                isnum = true
                for c_idx in children_idx
                    ix = children_arr[c_idx]
                    chdset = numvalued[ix]
                    isset &= chdset
                    if chdset
                        f_input[r] = setstorage[ix]
                    else
                        f_input[r] = numberstorage[ix]
                    end
                    r += 1
                end
                fval = MOI.eval_objective(evaluator, f_input)
                if isnum
                    numberstorage[k] = fval
                    println("AT k: $k,   START CALL assigment yeilds setstorage: $(numberstorage[k])")
                else
                    setstorage[k] = SetValuePost(x_values, fval, current_node)
                    println("AT k: $k,   START CALL assigment yeilds setstorage: $(setstorage[k])")
                end
                numvalued[k] = isnum
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == JuMP.Derivatives.CALLUNIVAR # univariate function
            op = nod.index
            println("univariate op: $op")
            child_idx = children_arr[adj.colptr[k]]
            chdset = numvalued[child_idx]
            if chdset
                child_val = numberstorage[child_idx]
            else
                child_val = setstorage[child_idx]
            end
            if op >= JuMP.Derivatives.USER_UNIVAR_OPERATOR_ID_START
                userop = op - JuMP.Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
                f = user_operators.univariate_operator_f[userop]
                fval = f(child_val)
            else
                #fval = JuMP.Derivatives.eval_univariate(op, child_val)
                fval = JuMP.Derivatives.eval_univariate(op, child_val)
            end
            if chdset
                @inbounds numberstorage[k] = fval
            else
                @inbounds setstorage[k] = SetValuePost(x_values, fval, current_node)
            end
            numvalued[k] = chdset
        elseif nod.nodetype == JuMP.Derivatives.COMPARISON
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            result = true
            for r in 1:n_children-1
                ix1 = children_arr[children_idx[r]]
                ix2 = children_arr[children_idx[r+1]]
                isnum1 = numvalued[ix1]
                isnum2 = numvalued[ix2]
                if isnum1
                    cval_lhs = numberstorage[ix1]
                else
                    cval_lhs = setstorage[ix1]
                end
                if isnum2
                    cval_rhs = numberstorage[ix2]
                else
                    cval_rhs = setstorage[ix2]
                end
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
            numberstorage[k] = result
        elseif nod.nodetype == JuMP.Derivatives.LOGIC
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            ix1 = children_arr[first(children_idx)]
            ix2 = children_arr[last(children_idx)]
            isnum1 = numvalued[ix1]
            isnum2 = numvalued[ix2]
            if isnum1
                cval_lhs = (numberstorage[ix1] == 1)
            else
                cval_lhs = (setstorage[ix1] == 1)
            end
            if isnum2
                cval_rhs = (numberstorage[ix2] == 1)
            else
                cval_rhs = (setstorage[ix2] == 1)
            end
            if op == 1
                numberstorage[k] = cval_lhs && cval_rhs
            elseif op == 2
                numberstorage[k] = cval_lhs || cval_rhs
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
        subexpr_values[k] = forward_eval(ex.setstorage, ex.numberstorage, ex.setvalued,
                                         ex.nd, ex.adj, ex.const_values,
                                         d.parameter_values, d.current_node,
                                         x, subexpr_values, user_input_buffer,
                                         user_operators=user_operators)
    end
    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.setvalued,
                     ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x, subexpr_values, user_input_buffer,
                     user_operators=user_operators)
    end
    for ex in d.constraints
        forward_eval(ex.setstorage, ex.numberstorage, ex.setvalued,
                     ex.nd, ex.adj, ex.const_values,
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
