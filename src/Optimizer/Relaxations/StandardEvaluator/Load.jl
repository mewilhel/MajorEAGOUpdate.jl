function copy_to_function(T::S,x::JuMP.FunctionStorage) where {S<:DataType}
    temp_set = Array{T}(undef,length(x.nd))
    temp_flt = Array{Float64}(undef,length(x.nd))
    temp_bool = Array{Bool}(undef,length(x.nd))
    FunctionSetStorage{T}(x.nd, x.adj, x.const_values, temp_set, temp_flt, temp_bool,
                          x.grad_sparsity, x.hess_I, x.hess_J, x.dependent_subexpressions)
end

function copy_to_subexpr(T::S,x::JuMP.SubexpressionStorage) where {S<:DataType}
    temp_set = Array{T}(undef,length(x.nd))
    temp_flt = Array{Float64}(undef,length(x.nd))
    temp_bool = Array{Bool}(undef,length(x.nd))
    SubexpressionSetStorage{T}(x.nd, x.adj, x.const_values, temp_set, temp_flt,
                               temp_bool, x.linearity)
end

function MinusObjective!(d::Evaluator{T}) where T<:Real
    if (d.has_nlobj)
        # shifts the adjacency matrix to introduce -f(x) as first element of nd array
        rowval = rowvals(d.objective.adj) .+ 1
        nzval = nonzeros(d.objective.adj)
        m, n = size(d.objective.adj)
        pushfirst!(d.objective.adj.colptr,1)
        d.objective.adj = SparseMatrixCSC{Bool,Int}(m+1,n+1,d.objective.adj.colptr,rowval,nzval)
        d.objective.adj[1,2] = true

        # shifts the node list (and parents)
        shift_nd = [JuMP.NodeData(JuMP.CALLUNIVAR,2,2)]
        for nd in d.objective.nd
            push!(shift_nd,JuMP.NodeData(nd.nodetype,nd.index,nd.parent+1))
        end
        d.objective.nd = shift_nd

        nvflag = length(d.objective.numvalued) > 0 ? d.objective.numvalued[1] : false
        pushfirst!(d.objective.numvalued,nvflag)
        pushfirst!(d.objective.numberstorage,0.0)
        pushfirst!(d.objective.setstorage,zero(T))
    end
end

# TO DO
function Build_NLP_Evaluator(S::R,src::T,x::Optimizer) where {R<:Type, T<:MOI.AbstractNLPEvaluator}

    # Checks to see if nlp data block evaluator is present
    if (typeof(src) != EAGO.EmptyNLPEvaluator)

        # Creates the empty evaluator
        d = Evaluator{S}(src.m)

        num_variables_ = JuMP.num_variables(d.m)
        d.variable_number = num_variables_
        nldata::JuMP.NLPData = d.m.nlp_data

        # Copy state of user-defined multivariable functions
        d.has_user_mv_operator = src.disable_2ndorder
        d.parameter_values = nldata.nlparamvalues
        d.last_x = fill(NaN, d.variable_number)

        # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
        d.disable_1storder = !in(:cv_grad,fieldnames(eltype(d)))
        d.disable_2ndorder = !in(:cv_hess,fieldnames(eltype(d)))

        # Add objective functions, constraints, subexpressions
        d.has_nlobj = isa(nldata.nlobj, JuMP.NonlinearExprData)
        if (src.has_nlobj)
            d.objective = copy_to_function(S,src.objective)
        end

        for i in 1:length(src.constraints)
            push!(d.constraints,copy_to_function(S,src.constraints[i]))
        end

        d.subexpression_order = src.subexpression_order
        d.subexpression_linearity = src.subexpression_linearity
        d.subexpressions_as_julia_expressions = Any[]
        if isdefined(src,:subexpressions_as_julia_expressions)
            d.subexpressions_as_julia_expressions = src.subexpressions_as_julia_expressions
        end

        d.subexpression_values_set = S[]
        d.subexpression_values_flt = Float64[]
        d.subexpressions = SubexpressionSetStorage{S}[]
        for i in 1:length(src.subexpressions)
            temp = copy_to_subexpr(S,src.subexpressions[i])
            push!(d.subexpressions,temp)
        end
        d.subexpression_values_set = fill(NaN,length(d.subexpressions))
        d.subexpression_values_flt = fill(NaN,length(d.subexpressions))

        # USER OUTPUT BUFFERS??????
        d.fw_repeats = x.EvalWalkRepts
        d.has_reverse = x.EvalReverse
        d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension)) # DO I NEED THIS

        return d
    end
end
