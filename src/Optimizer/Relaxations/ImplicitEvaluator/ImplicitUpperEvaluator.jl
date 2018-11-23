const RD_COMPILE_SWITCH = 3000

# FIX ME
mutable struct ImplicitUpperEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool

    fg::Function
    func_eval::Bool

    np::Int
    ny::Int
    ng::Int
    last_y::Vector{Float64}
    diff_result
    diff_tape

    value_storage::Vector{Float64}
    jacobian_storage::VecOrMat{Float64}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
    function ImplicitUpperEvaluator()
        d = new()
        d.disable_1storder = false
        d.has_nlobj = false
        d.func_eval = false
        return d
    end
end

function reform_1!(out,y,f,g,h,np,ng,nx)
    out[1] = f(y[1:np],y[(np+1):(np+nx)])
    out[2:(ng+1)] = g(y[1:np],y[(np+1):(np+nx)])
    out[(ng+2):(nx+ng+1)] = h(y[1:np],y[(np+1):(np+nx)])
    out[(2+ng+nx):(1+ng+2*nx)] = -out[(ng+2):(nx+ng+1)]
end

function reform_2!(out,y,g,h,np,ng,nx)
    out[1:ng] = g(y[1:np],y[(np+1):(np+nx)])
    out[(ng+1):(nx+ng)] = h(y[1:np],y[(np+1):(np+nx)])
    out[(1+ng+nx):(ng+2*nx)] = -h(y[1:np],y[(np+1):(np+nx)])
end

function reform_3!(out,y,f,h,np,nx)
    out[1] = f(y[1:np],y[(np+1):(np+nx)])
    out[2:(nx+1)] = h(y[1:np],y[(np+1):(np+nx)])
    out[(nx+2):(2*nx+1)] = -out[2:(nx+1)]
end

# LOOKS GREAT!
function build_upper_evaluator!(d::ImplicitUpperEvaluator; obj = nothing, constr = nothing,
                                impfun = nothing, np::Int = 0, nx::Int = 0, ng::Int = 0,
                                user_sparse::Vector{Tuple{Int64,Int64}} = nothing)
    d.ny = np + nx
    d.ng = ng
    if (d.ng > 0 && obj != nothing)
        d.fg = (out,x) -> reform_1!(out,x,obj,constr,impfun,np,ng,nx)
        d.has_nlobj = true
    elseif (obj == nothing)
        d.fg = (out,x) -> reform_2!(out,x,constr,impfun,np,ng,nx)
    else
        d.fg = (out,x) -> reform_3!(out,x,obj,impfun,np,nx)
        d.has_nlobj = true
    end
    d.value_storage = (obj != nothing) ? zeros(1+ng+2*nx) : zeros(ng+2*nx)
    d.diff_result = (obj != nothing) ? zeros(1+ng+2*nx,np+nx) : zeros(ng+2*nx,np+nx)
    d.last_y = zeros(d.ny)
    d.diff_tape = ReverseDiff.JacobianTape(d.fg, zeros(1+ng+2*nx), d.last_y)
    if length(d.diff_tape) < RD_COMPILE_SWITCH
        d.diff_tape = ReverseDiff.compile(d.diff_tape)
    end
end

# LOOKS GREAT!
function calc_functions!(d::ImplicitUpperEvaluator,y)
    if (d.last_y != y)
        if ~d.disable_1storder
            d.fg(d.value_storage,y)
            out = ReverseDiff.jacobian!(d.diff_tape, y)
            ReverseDiff.jacobian!(d.diff_result, d.diff_tape, y)
        end
        d.fg(d.value_storage,y)
        d.func_eval = true
    end
end

# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitUpperEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = 0.0
        if (d.has_nlobj)
            calc_functions!(d,y)
            val = d.value_storage[1]
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitUpperEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if (d.ng+d.nx) > 0
            calc_functions!(d,y)
            g[:] = d.value_storage[2:end]
        end
    end
    return
end
#=
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
=#
