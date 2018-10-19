MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.ZeroOne}) = true

MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true

MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

#=
    Add unconstrained variables to models
=#
function MOI.add_variable(m::Optimizer)
    m.VariableNumber += 1
    m.NonlinearVariable[m.VariableNumber] = false
    push!(m.VariableInfo, VariableInfo())
    return MOI.VariableIndex(length(m.VariableInfo))
end
MOI.add_variables(m::Optimizer, n::Int) = [MOI.add_variable(m) for i in 1:n]


function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, zo::MOI.ZeroOne)
    vi = v.variable
    check_inbounds(m, vi)
    has_upper_bound(m, vi) && error("Upper bound on variable $vi already exists.")
    has_lower_bound(m, vi) && error("Lower bound on variable $vi already exists.")
    is_fixed(m, vi) && error("Variable $vi is fixed. Cannot also set upper bound.")
    m.VariableInfo[vi.value].lower_bound = 0.0
    m.VariableInfo[vi.value].upper_bound = 1.0
    m.VariableInfo[vi.value].has_lower_bound = true
    m.VariableInfo[vi.value].has_upper_bound = true
    m.VariableInfo[vi.value].is_integer = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}(vi.value)
end

#=
    Add single variable constraints
=#
function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, lt::MOI.LessThan{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set upper bound.")
    end
    m.VariableInfo[vi.value].upper_bound = lt.upper
    m.VariableInfo[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, gt::MOI.GreaterThan{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set lower bound.")
    end
    m.VariableInfo[vi.value].lower_bound = gt.lower
    m.VariableInfo[vi.value].has_lower_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, eq::MOI.EqualTo{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Variable $vi has a lower bound. Cannot be fixed.")
    end
    if has_upper_bound(m, vi)
        error("Variable $vi has an upper bound. Cannot be fixed.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is already fixed.")
    end
    m.VariableInfo[vi.value].lower_bound = eq.value
    m.VariableInfo[vi.value].upper_bound = eq.value
    m.VariableInfo[vi.value].has_lower_bound = true
    m.VariableInfo[vi.value].has_upper_bound = true
    m.VariableInfo[vi.value].is_fixed = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, eq::MOI.Interval{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(eq.lower)
        error("Invalid fixed value $(gt.lower).")
    end
    if isnan(eq.upper)
        error("Invalid fixed value $(gt.upper).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set interval bounds.")
    end
    m.VariableInfo[vi.value].lower_bound = eq.lower
    m.VariableInfo[vi.value].upper_bound = eq.upper
    m.VariableInfo[vi.value].has_lower_bound = true
    m.VariableInfo[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(vi.value)
end

#=
    Linear and quadratic constraints
=#

macro define_addconstraint_linear(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(m, func)
            push!(m.$(array_name), (func, set, length(func.terms)))
            indx = MOI.ConstraintIndex{$function_type, $set_type}(length(m.$(array_name)))
            m.ConstraintConvexity[indx] = true
            return indx
        end
    end
end

macro define_addconstraint_quadratic(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(m, func)
            for i in func.affine_terms m.NonlinearVariable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.NonlinearVariable[i.variable_index_1.value] = true
                m.NonlinearVariable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set, length(m.$(array_name))+1))
            indx = MOI.ConstraintIndex{$function_type, $set_type}(length(m.$(array_name)))
            m.ConstraintConvexity[indx] = false
            return indx
        end
    end
end

@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.LessThan{Float64} LinearLEQConstraints
@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.GreaterThan{Float64} LinearGEQConstraints
@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.EqualTo{Float64} LinearEQConstraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.LessThan{Float64} QuadraticLEQConstraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.GreaterThan{Float64} QuadraticGEQConstraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.EqualTo{Float64} QuadraticEQConstraints

#=
    Adds nonlinear constraints
=#
function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    m.NLPData = nlp_data
    return
end
