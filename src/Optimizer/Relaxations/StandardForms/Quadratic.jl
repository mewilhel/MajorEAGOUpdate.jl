function Classify_Quadratic()
end

function Relax_Convex_Quadratic_Inner!(trg::T, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                 lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64}) where {T<:MOI.AbstractOptimizer}

    VarNum = length(n.LowerVar)
    temp_sto = zeros(Float64,VarNum)
    saf = ScalarAffineFunction{Float64}(ScalarAffineTerm{Float64}[],func.constant)
    for term in func.quadratic_terms
        @inbounds temp_constant = x0[VItoSto[term.variable_index_1]]
        @inbounds temp_constant *= x0[VItoSto[term.variable_index_2]]
        saf.constant -= term.coefficient*temp_constant
    end

    # Adds affine terms
    for term in func.affine_terms
        @inbounds temp_sto[VItoSto[term.variable_index]] += term.coefficient
    end

    # Adds quadratic terms
    for term in func.quadratic_terms
        @inbounds temp_sto[VItoSto[term.variable_index_1]] += 2.0*term.coefficient*x0[VItoSto[term.variable_index_1]]
    end

    # Adds nonzero terms
    temp = 0.0
    for i=1:VarNum
        @inbounds temp = temp_sto[i]
        if (temp != 0.0)
            push!(saf.terms,ScalarAffineTerm{Float64}(temp,StoToVi[i]))
        end
    end
    if (lower == upper)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
        saf.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    elseif (lower != -Inf)
        saf.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    else
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
    end
end

function Relax_Nonconvex_Quadratic!(trg::T, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                    lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64}) where {T<:MOI.AbstractOptimizer}
    saf = ScalarAffineFunction{Float64}(ScalarAffineTerm{Float64}[],func.constant)
    VarNum = length(n.LowerVar)
    quadratic_coefficient = zeros(Float64,VarNum)
    quadratic_constant = func.constant
    for term in func.quadratic_terms
        a = term.coefficient
        @inbounds xL1 = n.LowerVar[VItoSto[term.variable_index_1]]
        @inbounds xU1 = n.UpperVar[VItoSto[term.variable_index_1]]
        if VItoSto[term.variable_index_1] == VItoSto[term.variable_index_2]               # quadratic terms
            if (a > 0.0)
                @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] = a*(2.0*xL1+1.0)
                @inbounds quadratic_constant -= a*xL1*xU1
            else
                @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] = a*(xL1 + xU1)
                @inbounds quadratic_constant -= a*xL1*xU1
            end
        else
            @inbounds xL2 = n.LowerVar[VItoSto[term.variable_index_2]]
            @inbounds xU2 = n.UpperVar[VItoSto[term.variable_index_2]]
            if (a > 0.0)
                @inbounds check = (xU1 - xL1)*x0[VItoSto[term.variable_index_2]] + (xU2 - xL2)*x0[VItoSto[term.variable_index_1]] - xU1*xU2 + xL1*xL2 <= 0.0
                if check
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] += a*xL2
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_2]] += a*xL1
                    @inbounds quadratic_constant -= a*xL1*xL2
                else
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] += a*xU2
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_2]] += a*xU1
                    @inbounds quadratic_constant -= a*xU1*xU2
                end
            else
                @inbounds check = (xU1 - xL1)*x0[VItoSto[term.variable_index_2]] + (xU2 - xL2)*x0[VItoSto[term.variable_index_1]] - xU1*xL2 + xL1*xU2 <= 0.0
                if check
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] += a*xL2
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_2]] += a*xU1
                    @inbounds quadratic_constant -= a*xU1*xL2
                else
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_1]] += a*xU2
                    @inbounds quadratic_coefficient[VItoSto[term.variable_index_2]] += a*xL1
                    @inbounds quadratic_constant -= a*xL1*xU2
                end
            end
        end
    end

    temp = 0.0
    for i=1:VarNum
        @inbounds temp = quadratic_coefficient[i]
        if (temp != 0.0)
            push!(saf.terms,ScalarAffineTerm{Float64}(temp,StoToVi[i]))
        end
    end

    if (lower == upper)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
        saf.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        saf2.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    elseif (lower != -Inf)
        saf.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        saf.constant *= -1.0
        for i in 1:length(saf.terms)
            @inbounds saf.terms[i] *= -1.0
        end
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    else
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
    end
end

function RelaxQuadratic!(trg::T, src::Optimizer, n::NodeBB) where {T<:MOI.AbstractOptimizer}
    # Relax Convex Quadratic Terms
    for (func, set, ind) in src.QuadraticLEQConstraints
        if QuadraticConvexity[ind]
            Relax_Convex_Quadratic!(trg, src, func, -Inf, set.upper, n)
        else
            Relax_Nonconvex_Quadratic!(trg, src, func, -Inf, set.upper, n)
        end
    end
    for (func,set,ind) in src.QuadraticGEQConstraints
        if QuadraticConvexity[ind]
            Relax_Convex_Quadratic!(trg, src, func, set.lower, Inf, n)
        else
            Relax_Nonconvex_Quadratic!(trg, src, func, set.lower, Inf, n)
        end
    end
    for (func,set,ind) in src.QuadraticEQConstraints
        if QuadraticConvexity[ind]
            Relax_Convex_Quadratic!(trg, src, func, set.value, set.value, n)
        else
            Relax_Nonconvex_Quadratic!(trg, src, func, set.value, set.value, n)
        end
    end
end
