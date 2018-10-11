function convert(FunctionSetStorage{T},x::FunctionStorage) where T
    FunctionSetStorage(x.nd, x.adj, x.const_values, T[], x.grad_sparsity, x.hess_I, x.hess_J)
end

function convert(SubexpressionSetStorage{T},x::FunctionStorage) where T
    SubexpressionSetStorage(x.nd, x.adj, x.const_values, T[], x.Linearity)
end

# TO DO
function Build_NLP_Evaluator(src<:AbstractNLPEvaluator)

    # Creates the empty evaluator
    d =  Evaluator()

    num_variables_ = num_variables(d.m)
    nldata::NLPData = d.m.nlp_data

    # Copy state of user-defined multivariable functions
    d.has_user_mv_operator = src.disable_2ndorder
    d.parameter_values = nldata.nlparamvalues
    d.last_x = fill(NaN, d.variable_number)

    d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension)) # DO I NEED THIS

    # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
    d.disable_1storder = !in(:cv_grad,fieldnames(d))
    d.disable_2ndorder = !in(:cv_hess,fieldnames(d))
    for feat in requested_features
        if !(feat in MOI.features_available(d))
            error("Unsupported feature $feat")
        end
    end

    # Add objective functions, constraints, subexpressions
    d.has_nlobj = isa(nldata.nlobj, NonlinearExprData)
    d.objective = convert(FunctionSetStorage,src.objective)
    for i in 1:length(src.constraints)
        push!(d.constraints,convert(FunctionSetStorage,src.constraints[i]))
    end
    for i in 1:length(src.subexpressions)
        push!(d.subexpressions,convert(SubexpressionSetStorage,src.subexpressions[i]))
    end

    # copy the hessian sparsity pattern
    if !d.disable_2ndorder
        d.hessian_sparsity = src.hessian_sparsity
    end

    return d
end
