```
```
#=
TO DO: ADD CHEAP IPOPT SOLVE THAT DOESN'T REQUIRE REBUILDING MOI MODEL
struct IpoptStorage <: SubProblemInfo
    num_variables
    num_linear_le_constraints
    num_linear_ge_constraints
    num_linear_eq_constraints
    nlp_row_offset
    num_quadratic_constraints
    num_nlp_constraints
    num_constraints
end


function CreateSolverInternal(model::Optimizer,::Ipopt.Optimizer)
    IpoptStorage.num_variables = length(model.variable_info)
    IpoptStorage.num_linear_le_constraints = length(model.linear_le_constraints)
    IpoptStorage.num_linear_ge_constraints = length(model.linear_ge_constraints)
    IpoptStorage.num_linear_eq_constraints = length(model.linear_eq_constraints)
    IpoptStorage.nlp_row_offset = nlp_constraint_offset(model)
    IpoptStorage.num_quadratic_constraints = nlp_constraint_offset(model) -
                                quadratic_le_offset(model)
    IpoptStorage.num_nlp_constraints = length(model.nlp_data.constraint_bounds)
    IpoptStorage.num_constraints = IpoptStorage.num_nlp_constraints + IpoptStorage.nlp_row_offset

    IpoptStorage.evaluator = model.nlp_data.evaluator
    IpoptStorage.features = MOI.features_available(evaluator)
    IpoptStorage.has_hessian = (:Hess in features)
    IpoptStorage.init_feat = [:Grad]
    IpoptStorage.has_hessian && push!(init_feat, :Hess)
    IpoptStorage.num_nlp_constraints > 0 && push!(init_feat, :Jac)

    MOI.initialize(evaluator, init_feat)
    jacobian_sparsity = jacobian_structure(model)
    hessian_sparsity = has_hessian ? hessian_lagrangian_structure(model) : []

    # Objective callback
    if model.sense == MOI.MinSense
        objective_scale = 1.0
    elseif model.sense == MOI.MaxSense
        objective_scale = -1.0
    else # FeasibilitySense
        # TODO: This could produce confusing solver output if a nonzero
        # objective is set.
        objective_scale = 0.0
    end

    eval_f_cb(x) = objective_scale * eval_objective(model, x)

    # Objective gradient callback
    function eval_grad_f_cb(x, grad_f)
        eval_objective_gradient(model, grad_f, x)
        Compat.rmul!(grad_f,objective_scale)
    end

    # Constraint value callback
    eval_g_cb(x, g) = eval_constraint(model, g, x)

    # Jacobian callback
    function eval_jac_g_cb(x, mode, rows, cols, values)
        if mode == :Structure
            for i in 1:length(jacobian_sparsity)
                rows[i] = jacobian_sparsity[i][1]
                cols[i] = jacobian_sparsity[i][2]
            end
        else
            eval_constraint_jacobian(model, values, x)
        end
    end

    if has_hessian
        # Hessian callback
        function eval_h_cb(x, mode, rows, cols, obj_factor,
            lambda, values)
            if mode == :Structure
                for i in 1:length(hessian_sparsity)
                    rows[i] = hessian_sparsity[i][1]
                    cols[i] = hessian_sparsity[i][2]
                end
            else
                obj_factor *= objective_scale
                eval_hessian_lagrangian(model, values, x, objective_scale*obj_factor, lambda)
            end
        end
    else
        eval_h_cb = nothing
    end

    x_l = [v.lower_bound for v in model.variable_info]
    x_u = [v.upper_bound for v in model.variable_info]

    constraint_lb, constraint_ub = constraint_bounds(model)

    model.inner.x = [v.start for v in model.variable_info]
end

function MOI.optimize!(model::Optimizer,::Ipopt.Optimizer)

    model.inner = createProblem(num_variables, x_l, x_u, num_constraints,
                            constraint_lb, constraint_ub,
                            length(jacobian_sparsity),
                            length(hessian_sparsity),
                            eval_f_cb, eval_g_cb, eval_grad_f_cb, eval_jac_g_cb,
                            eval_h_cb)

.
    addOption(model.inner, "check_derivatives_for_naninf", "yes")

    if !IpoptStorage.has_hessian
        addOption(model.inner, "hessian_approximation", "limited-memory")
    end
    if num_nlp_constraints == 0 && num_quadratic_constraints == 0
        addOption(model.inner, "jac_c_constant", "yes")
        addOption(model.inner, "jac_d_constant", "yes")
        if !model.nlp_data.has_objective
            addOption(model.inner, "hessian_constant", "yes")
        end
    end

    for (name,value) in model.options
        sname = string(name)
        if match(r"(^resto_)", sname) != nothing
            sname = replace(sname, r"(^resto_)", "resto.")
        end
        addOption(model.inner, sname, value)
    end
    solveProblem(model.inner)
end
=#
