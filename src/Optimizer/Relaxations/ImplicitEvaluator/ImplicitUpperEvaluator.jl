# state gradient  functions and expand forwarddiff to allow for computation of state then f,g
mutable struct ImplicitUpperEvaluatorEvaluator{T<:Real} <: MOI.AbstractNLPEvaluator
    current_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool

    var_number::Int
    cnstr_number::Int
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
    function Evaluator{T}(m) where T<:Real
        d = new()
        return d
    end
end

function SetupEvaluator(d,f,g)
    if (cnstr_number > 0) && has_nlobj
        d.fg = x-> vcat(f(x),g(x))
    elseif has_nlobj
        d.fg = x-> f(x)
    else
        d.fg = x-> g(x)
    end
    d.diff_result = DiffResults.JacobianResult(x)
end


# LOOKS GREAT!
function calc_functions!(d::ImplicitUpperEvaluator,p)
    if (d.last_p != p)
        if ~disable_1storder
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
        val = zero(eltype(p))
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
        if d.num_constraints > 0
            calc_functions!(d,p)
            for i in 1:d.num_constraints
                g[i] = cnstr_relax[i].cv
            end
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitUpperEvaluator, df, p)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            calc_functions!(d,p)
            for j in 1:d.num_control
                df[j] = d.obj_relax.cv_grad[j]
            end
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
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.num_constraints for idx in 1:d.num_controls]
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

# TO DO
function MOI.eval_constraint_jacobian(d::ImplicitUpperEvaluator,g,p)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        forward_reverse_pass(d,p)
        #t = typeof(d.constraints[1].setstorage[1])
        #g = zero.(g)
        for i in 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j in 1:length(d.constraints[i].setstorage[1].cv_grad)
                    g[i,j] = d.constraints[i].setstorage[1].cv_grad[j]
                end
            else
                for j in 1:length(d.constraints[i].setstorage[1].cv_grad)
                    g[i,j] = 0.0
                end
            end
        end
    #end
    return
end

# TO DO
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

# TO DO
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
