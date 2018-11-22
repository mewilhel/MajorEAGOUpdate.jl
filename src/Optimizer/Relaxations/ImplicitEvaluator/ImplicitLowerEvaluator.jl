mutable struct ImplicitLowerEvaluator{N} <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
    last_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    objective::Function
    np::Int
    nx::Int
    ng::Int
    constraints::Function
    jacobian_sparsity::Vector{Tuple{Int64,Int64}}
    obj_eval::Bool
    cnstr_eval::Bool

    last_p::Vector{Float64}
    ref_p::Vector{Float64}
    var_relax::Vector{MC{N}}
    state_relax::Vector{MC{N}}
    state_ref_relaxation::Vector{Vector{MC{N}}}
    obj_relax::MC{N}
    cnstr_relax::Vector{MC{N}}
    imp_opts::mc_opts

    # timer
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function ImplicitLowerEvaluator{N}() where N
        d = new()

        d.objective = nothing
        d.constraints = nothing
        d.obj_eval = false
        d.cnstr_eval = false
        d.jacobian_sparsity = nothing

        d.np = 0
        d.nx = 0
        d.ng = 0

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        d.var_relax = [zero(MC{N}) for i in 1:N]
        d.state_relax = [zero(MC{N})]
        d.state_ref_relaxation = [[zero(MC{N})]]
        d.obj_relax = zero(MC{N})
        d.cnstr_relax = [zero(MC{N})]

        last_p = Float64[0.0]
        ref_p = Float64[0.0]

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false

        d.current_node = NodeBB()
        d.last_node = NodeBB()

        d.imp_opts = mc_opts()

        return d
    end
end

function build_lower_evaluator!(d::ImplicitLowerEvaluator; obj = nothing, constr = nothing,
                                impfun = nothing, np::Int = 0, nx::Int = 0, ng::Int = 0,
                                user_sparse::Vector{Tuple{Int64,Int64}} = nothing)

    # setup objective and constraint functions
    d.objective = obj
    d.constraints = constr

    # has a nonlinear objective?
    (obj == nothing) || (d.has_nlobj = true)

    # get dimension sizes
    d.np = np; d.nx = nx;
    d.ng = (constr == nothing) ? 0 : length(constr(zeros(nx),zeros(np)))

    # set implicit routine information
    d.imp_opts.nx = nx
    d.imp_opts.np = np

    # preallocates the storage variables
    d.state_relax = MC{N}[zero(MC{N}) for i in 1:nx]
    d.cnstr_relax = MC{N}[zero(MC{N}) for i in 1:ng]
    d.state_ref_relaxation = Vector{MC{N}}[MC{N}[zero(MC{N}) for i in 1:nx] for j in 1:d.imp_opts.kmax]

    # allocates the reference points
    d.last_p = zeros(Float64,np)
    d.ref_p = zeros(Float64,np)

    if (user_sparse == nothing)
        sparse_pattern = Tuple{Int64,Int64}[]
        for i in 1:nx
            for j in 1:np
                push!(sparse_pattern,(i,j))
            end
        end
        d.jacobian_sparsity = sparse_pattern
    end
end

# LOOKS GREAT!
function relax_implicit!(d::ImplicitLowerEvaluator,p)
    # Generate new parameters for implicit relaxation if necessary
    if d.current_node != d.last_node
        d.obj_eval = false
        d.cnstr_eval = false
        for i in 1:np
            d.ref_p[i] = (d.current_node.LowerVar[i]+d.current_node.UpperVar[i])/2.0
            d.P[i] = IntervalType(d.current_node.LowerVar[i],d.current_node.UpperVar[i])
        end
        for j in (np+1):(np+nx)
            shiftj = j - np
            d.X[shiftj] = IntervalType(d.current_node.LowerVar[j],d.current_node.UpperVar[j])
        end
        d.state_ref_relaxation = GenExpansionParams(d.h, d.hj, d.X, d.P, d.ref_p, d.mc_opts)
    end
    # Generate new value of implicit relaxation
    if d.ref_p != p
        d.obj_eval = false
        d.cnstr_eval = false
        pMC = MC{np}(p,d.P)
        d.state_relax = MC_impRelax(d.h, d.hj, pMC, d.ref_p, d.X, d.P, d.mc_opts, d.state_ref_relaxation)
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
            val = d.obj_relax.cv
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitLowerEvaluator, g, p)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,p)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.cnstr_relax[i].cv
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
            for j in 1:d.np
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
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.num_constraints for idx in 1:d.np]
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

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitLowerEvaluator,g,p)
    #d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            relax_implicit!(d,p)
            relax_constraints!(d)
            g = zero.(g)
            for (i,j) in d.jacobian_sparsity
                g[i,j] = d.cnstr_relax[i].cv_grad[j]
            end
        end
    #end
    return
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian_product(d::ImplicitLowerEvaluator, y, p, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_implicit!(d,p)
                relax_constraints!(d)
                y = zero.(y)
                for (i,j) in d.jacobian_sparsity
                    y[i] += d.cnstr_relax[i].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# PROBABLY DONE
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitLowerEvaluator, y, p, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_implicit!(d,p)
                relax_constraints!(d)
                y = zero.(y)
                for (i,j) in d.jacobian_sparsity
                    y[i] += d.cnstr_relax[i].cv_grad[j]*w[j]
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
