# FIX ME
mutable struct ImplicitUpperEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool

    np::Int
    ng::Int
    last_p::Vector{Float64}
    diff_result::Union{DiffResults.JacobianResult,DiffResults.HessianResult}
    diff_config::Union{ForwardDiff.JacobianConfig,ForwardDiff.HessianConfig}

    value_storage::Vector{Float64}
    jacobian_storage::VecOrMat{Float64}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
    function ImplicitUpperEvaluator()
        d = new()
        return d
    end
end

# LOOKS GREAT!
function SetupEvaluator(d,f,g)
    if (d.ng > 0) && d.has_nlobj
        d.fg = x -> vcat(f(x),g(x))
    elseif has_nlobj
        d.fg = x -> f(x)
    else
        d.fg = x -> g(x)
    end
    d.diff_result = DiffResults.JacobianResult(x)
end


# LOOKS GREAT!
function calc_functions!(d::ImplicitUpperEvaluator,p)
    if (d.last_p != p)
        if ~d.disable_1storder
            ForwardDiff.jacobian!(d.diff_result, d.fg, p, d.diff_config)
            d.value_storage = DiffResults.value(d.diff_result)
            d.jacobian_storage = DiffResults.jacobian(d.diff_result)
        else
            d.value_storage = d.fg(p)
        end
        d.func_eval = true
    end
end

# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitUpperEvaluator, p)
    d.eval_objective_timer += @elapsed begin
        val = 0.0
        if (d.has_nlobj)
            calc_functions!(d,p)
            val = d.value_storage[1]
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitUpperEvaluator, g, p)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            calc_functions!(d,p)
            g[:] = d.value_storage[2:end]
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitUpperEvaluator, df, p)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            calc_functions!(d,p)
            df[:] = d.jacobian_storage[1,1:d.np]
        else
            error("No nonlinear objective.")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.jacobian_structure(d::ImplicitUpperEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else # else assume dense pattern
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.ng for idx in 1:d.np]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitUpperEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitUpperEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitUpperEvaluator, dg, p)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            calc_functions!(d,p)
            dg[:,:] = jacobian_storage[2:end,:]
        end
    #end
    return
end

# FIX ME
function MOI.eval_constraint_jacobian_product(d::ImplicitUpperEvaluator, y, p, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,p)
            t = typeof(d.constraints[1].setstorage[1])
            y = zeros(t,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
            for i in 1:length(d.constraints)
                if ~d.constraints[i].numvalued[1]
                    for j in 1:d.variable_number
                        y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                    end
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# FIX ME
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitUpperEvaluator, y, p, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            forward_reverse_pass(d,p)
            #t = typeof(d.constraints[1].setstorage[1])
            y = zeros(Float64,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
            for i in 1:length(d.constraints)
                if ~d.constraints[i].numvalued[1]
                    for j in 1:d.variable_number
                        y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                    end
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# LOOKS GREAT
function MOI.features_available(d::ImplicitUpperEvaluator)
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

# LOOKS GREAT
function MOI.initialize(d::ImplicitUpperEvaluator, requested_features::Vector{Symbol}) end

# LOOKS GREAT
MOI.objective_expr(d::ImplicitUpperEvaluator) = error("EAGO.ImplicitUpperEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitUpperEvaluator) = error("EAGO.ImplicitUpperEvaluator doesn't provide expression graphs of constraint functions.")
