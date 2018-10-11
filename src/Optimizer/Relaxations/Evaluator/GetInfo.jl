# looks good
function MOI.eval_objective(d::Evaluator, x)
    d.eval_objective_timer += @elapsed begin
        forward_reverse_pass(d,x)
        val = zero(eltype(x))
        if d.has_nlobj
            val = d.objective.storage[1].cv
        else
            error("No nonlinar objective.")
        end
    end
    return val
end

# looks good
function MOI.eval_constraint(d::Evaluator, g, x)
    d.eval_constraint_timer += @elapsed begin
        forward_reverse_pass(d,x)
        for i in 1:length(d.constraints)
            g[i] = d.constraints[i].forward_storage[1].cv
        end
    end
    return
end

# looks good
function MOI.eval_objective_gradient(d::Evaluator, df, x)
    d.eval_objective_timer += @elapsed begin
        forward_reverse_pass(d,x)
        val = zero(eltype(x))
        if d.has_nlobj
            val = d.objective.forward_storage[1].cv_grad
        else
            error("No nonlinar objective.")
        end
    end
    return val
end

# looks good
function MOI.jacobian_structure(d::Evaluator)
    jacobian_sparsity = Tuple{Int64,Int64}[]
    for row in 1:length(d.constraints)
        row_sparsity = d.constraints[row].grad_sparsity
        for idx in row_sparsity
            push!(jacobian_sparsity, (row, idx))
        end
    end
    return jacobian_sparsity
end

# looks good
function MOI.hessian_lagrangian_structure(d::Evaluator)
    !d.disable_2ndorder || error("Hessian computations were not requested on the call to initialize!.")
    return d.hessian_sparsity
end

# looks good
function _hessian_lagrangian_structure(d::Evaluator)
    hessian_sparsity = Tuple{Int64,Int64}[]
    if d.has_nlobj
        for idx in 1:length(d.objective.hess_I)
            push!(hessian_sparsity, (d.objective.hess_I[idx], d.objective.hess_J[idx]))
        end
    end
    for ex in d.constraints
        for idx in 1:length(ex.hess_I)
            push!(hessian_sparsity, (ex.hess_I[idx], ex.hess_J[idx]))
        end
    end
    return hessian_sparsity
end

# looks good
function MOI.eval_constraint_jacobian(d::Evaluator)
    d.eval_constraint_jacobian_timer += @elapsed begin
        forward_reverse_pass(d,x)
        t = typeof(d.constraints[i].forward_storage[1])
        g = zeros(t,length(d.constraints[i].forward_storage[1].cv_grad),length(d.constraints))
        for i in 1:length(d.constraints)
            g[:,i] = d.constraints[i].forward_storage[1].cv_grad
        end
    end
    return
end

# TO DO (CHECK GRADIENT DIMS)
function MOI.eval_constraint_jacobian_product(d::Evaluator, y, x, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            t = typeof(d.constraints[i].forward_storage[1])
            g = zeros(t,length(d.constraints[i].forward_storage[1].cv_grad),length(d.constraints))
            for i in 1:length(d.constraints)
                y[i] = d.constraints[i].forward_storage[1].cv_grad*w
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            t = typeof(d.constraints[i].forward_storage[1])
            y = zero(y)
            for i in 1:length(d.constraints)
                for j in 1:d.variable_number
                    y[i] += d.constraints[i].forward_storage[1].cv_grad[i]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_hessian_lagrangian_product(d::Evaluator, h, x, v, σ, μ)

    return h
end

# TO DO
function MOI.eval_hessian_lagrangian(d::Evaluator, H, x, σ, μ)
    if (!d.disable_2ndorder)
        d.eval_hessian_lagrangian_timer += @elapsed begin
            forward_reverse_pass(d,x)
            t = typeof(d.constraints[i].forward_storage[1])
            H[:,:] =  σ*d.objective.storage[1].cv_hess
            for i in 1:length(d.constraints)
                H[:,:] += μ[i]*d.constraints[i].forward_storage[1].cv_hess
            end
        end
    end
end

# looks good
function MOI.features_available(d::Evaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    if !d.disable_2ndorder
        push!(features,:Hess)
        push!(features,:HessVec)
    end
    return features
end

# looks good, doesn't do anything, EAGO builds the evaluator and attaches it to lower problems
function MOI.initialize(d::Evaluator, requested_features::Vector{Symbol}) end

# looks good
MOI.objective_expr(d::Evaluator) = error("EAGO.Evaluator doesn't provide expression graphs of constraint functions.")
# looks good
MOI.constraint_expr(d::Evaluator) = error("EAGO.Evaluator doesn't provide expression graphs of constraint functions.")
