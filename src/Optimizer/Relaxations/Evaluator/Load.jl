function convert(::FunctionSetStorage{T},x::JuMP.FunctionStorage) where T
    FunctionSetStorage(x.nd, x.adj, x.const_values, T[], x.grad_sparsity, x.hess_I, x.hess_J)
end

function convert(::SubexpressionSetStorage{T},x::JuMP.SubexpressionStorage) where T
    SubexpressionSetStorage(x.nd, x.adj, x.const_values, T[], x.Linearity)
end

# TO DO
function Build_NLP_Evaluator(S::R,src::T) where {R<:Type, T<:MOI.AbstractNLPEvaluator}

    # Creates the empty evaluator
    d = Evaluator{S}(src.m)

    num_variables_ = JuMP.num_variables(d.m)
    nldata::JuMP.NLPData = d.m.nlp_data

    # Copy state of user-defined multivariable functions
    d.has_user_mv_operator = src.disable_2ndorder
    d.parameter_values = nldata.nlparamvalues
    d.last_x = fill(NaN, d.variable_number)

    d.jac_storage = Array{Float64}(undef,max(num_variables_, d.m.nlp_data.largest_user_input_dimension)) # DO I NEED THIS

    # Set valued operators must contain a (sub)gradient and/or (sub)hessian field to support 1st/2nd order eval
    d.disable_1storder = !in(:cv_grad,fieldnames(eltype(d)))
    d.disable_2ndorder = !in(:cv_hess,fieldnames(eltype(d)))

    # Add objective functions, constraints, subexpressions
    d.has_nlobj = isa(nldata.nlobj, JuMP.NonlinearExprData)
    d.objective = convert(FunctionSetStorage,src.objective)
    #=
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
    =#
    return d
end
