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
                               temp_bool, x.Linearity)
end

# TO DO
function Build_NLP_Evaluator(S::R,src::T,x::Optimizer) where {R<:Type, T<:MOI.AbstractNLPEvaluator}

    # Checks to see if nlp data block evaluator is present
    if (typeof(src) != MajorEAGOUpdate.EmptyNLPEvaluator)

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
        c1 = d.has_nlobj
        c2 = d.objective
        c3 = src.objective

        d.has_nlobj = isa(nldata.nlobj, JuMP.NonlinearExprData)
        c1 = d.has_nlobj
        copied_func = copy_to_function(S,src.objective)
        d.objective = copy_to_function(S,src.objective)

        for i in 1:length(src.constraints)
            push!(d.constraints,copy_to_function(S,src.constraints[i]))
        end

        d.subexpression_order = src.subexpression_order
        d.subexpression_linearity = src.subexpression_linearity
        d.subexpressions_as_julia_expressions = Any[]
        if isdefined(src,:subexpressions_as_julia_expressions)
            d.subexpressions_as_julia_expressions = src.subexpressions_as_julia_expressions
        end

        d.subexpression_values = S[]
        d.subexpressions = SubexpressionSetStorage{S}[]
        for i in 1:length(src.subexpressions)
            push!(d.subexpressions,copy_to_subexpr(SubexpressionSetStorage{S},src.subexpressions[i]))
        end
        d.subexpression_values = fill(NaN,length(d.subexpressions))

        # USER OUTPUT BUFFERS??????
        d.fw_repeats = x.EvalWalkRepts
        d.has_reverse = x.EvalReverse
        d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension)) # DO I NEED THIS

        # copy the hessian sparsity pattern
        if !d.disable_2ndorder
            d.hessian_sparsity = src.hessian_sparsity
        end

        return d
    end
end
