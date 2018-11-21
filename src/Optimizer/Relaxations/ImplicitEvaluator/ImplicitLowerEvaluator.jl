mutable struct ImplicitLowerEvaluatorEvaluator{T<:Real} <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
    last_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    objective::Function
    num_control::Int
    num_constraints::Int
    constraints::Vector{Function}
    jacobian_sparsity::Vector{Tuple{Int64,Int64}}
    obj_eval::Bool
    cnstr_eval::Bool

    last_p::Vector{Float64}
    ref_p::Vector{Float64}
    var_relax::Vector{MC}
    state_relax::Vector{MC}
    state_ref_relaxation::Vector{Vector{MC}}
    obj_relax::MC
    cnstr_relax::Vector{MC}

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

# TO DO
function relax_implicit!(d::ImplicitLowerEvaluator,p)
    # Generate new parameters for implicit relaxation if necessary
    if current_node == last_node
        #GenExpansionParams(h::Function, hj::Function, X::Vector{IntervalType}, P::Vector{IntervalType}, pmid::Vector{Float64},mc_opts::mc_opts)
        d.state_ref_relaxation = GenExpansionParams(h, hj, X, P, pmid, mc_opts)
    end
    # Generate new value of implicit relaxation
    if d.ref_p != p
        d.ref_p = p
        d.obj_eval = false
        d.cnstr_eval = false
        # MC_impRelax(h::Function, hj::Function, p::Vector{MC{N}}, pmid::Vector{Float64}, X::Vector{IntervalType}, P::Vector{IntervalType}, mc_opts::mc_opts,param::Vector{Vector{MC{N}}}) where N
        d.state_relax = MC_impRelax(h, hj, p, pmid, X, P, mc_opts, d.state_ref_relaxation)
    else
        d.state_relax = d.state_ref_relaxation
    end
end

# LOOKS GREAT!
function relax_objective!(d::ImplicitLowerEvaluator)
    if ~d.obj_eval
        d.obj_relax = d.objective(d.state_relax,d.var_relax)
        d.obj_eval = true
    end
end

# LOOKS GREAT!
function relax_constraints!(d::ImplicitLowerEvaluator)
    if ~d.cntr_eval
        d.cnstr_relax = d.constraints(d.state_relax,d.var_relax)
        d.cnstr_eval = true
    end
end

# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitLowerEvaluator, p)
    d.eval_objective_timer += @elapsed begin
        val = zero(eltype(p))
        if d.has_nlobj
            relax_implicit!(d,p)
            relax_objective!(d)
            val = obj_relax.cv
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitLowerEvaluator, g, p)
    d.eval_constraint_timer += @elapsed begin
        if d.num_constraints > 0
            relax_implicit!(d,p)
            relax_constraints!(d)
            for i in 1:d.num_constraints
                g[i] = cnstr_relax[i].cv
            end
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitLowerEvaluator, df, p)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            relax_implicit!(d,p)
            relax_objective!(d)
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
function MOI.jacobian_structure(d::ImplicitLowerEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else # else assume dense pattern
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.num_constraints for idx in 1:d.num_controls]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitLowerEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitLowerEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# TO DO
function MOI.eval_constraint_jacobian(d::ImplicitLowerEvaluator,g,p)
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
function MOI.eval_constraint_jacobian_product(d::ImplicitLowerEvaluator, y, p, w)
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
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitLowerEvaluator, y, p, w)
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
function MOI.features_available(d::ImplicitLowerEvaluator)
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
function MOI.initialize(d::ImplicitLowerEvaluator, requested_features::Vector{Symbol}) end

# LOOKS GREAT
MOI.objective_expr(d::ImplicitLowerEvaluator) = error("EAGO.ImplicitLowerEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitLowerEvaluator) = error("EAGO.ImplicitLowerEvaluator doesn't provide expression graphs of constraint functions.")
