function Quadratic_Convexity!(src::Optimizer)
    for (func, set, ind) in src.QuadraticLEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(ind)
        src.ConstraintConvexity[MOIindx] = Is_Convex_Quadratic(func,src,1.0)
    end
    for (func,set,ind) in src.QuadraticGEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(ind)
        src.ConstraintConvexity[MOIindx] = Is_Convex_Quadratic(func,src,-1.0)
    end
    for (func,set,ind) in src.QuadraticEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(ind)
        src.ConstraintConvexity[MOIindx] = false
    end
end

function Is_Convex_Quadratic(func::MOI.ScalarQuadraticFunction{Float64},src::Optimizer,mult::Float64)
    NumVar = src.VariableNumber
    # Riguous Convexity Test
    Q = spzeros(NumVar,NumVar)
    for term in func.quadratic_terms
        if term.coefficient != 0.0
            Q[term.variable_index_1.value,term.variable_index_2.value] = mult*term.coefficient
        end
    end
    if length(Q.nzval) > 1
        eigmin = LinearAlgebra.eigmin(Array(Q))
        if (eigmin) > 0.0
            return true
        end
    else
        if Q.nzval[1] > 0.0
            return true
        else
            return false
        end
    end
end

function Relax_Convex_Quadratic_Inner!(trg, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                 lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64})

    VarNum = length(n.LowerVar)
    temp_sto = zeros(Float64,VarNum)
    terms_coeff = Float64[]
    terms_index = Int[]
    term_constant = func.constant
    for term in func.quadratic_terms
        temp_constant = x0[src.VItoSto[term.variable_index_1.value]]
        temp_constant *= x0[src.VItoSto[term.variable_index_2.value]]
        term_constant -= term.coefficient*temp_constant
    end

    # Adds affine terms
    for term in func.affine_terms
        temp_sto[src.VItoSto[term.variable_index.value]] += term.coefficient
    end

    # Adds quadratic terms
    for term in func.quadratic_terms
        @inbounds temp_sto[src.VItoSto[term.variable_index_1.value]] += 2.0*term.coefficient*x0[src.VItoSto[term.variable_index_1.value]]
    end

    temp = 0.0
    for (i,term) in enumerate(temp_sto)
        if (term != 0.0)
            push!(terms_coeff,term)
            push!(terms_index,src.StoToVI[i])
        end
    end

    varIndx = [MOI.VariableIndex(src.VItoSto[i]) for i in terms_index]

    if (lower == upper)
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
        term_constant *= -1.0
        for i in 1:length(terms_coeff)
            @inbounds terms_coeff[i] *= -1.0
        end
        saf1_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf1 = MOI.ScalarAffineFunction{Float64}(saf1_term,term_constant)
        MOI.add_constraint(trg,saf1,MOI.LessThan{Float64}(-lower))
    elseif (lower != -Inf)
        term_constant *= -1.0
        for i in 1:length(terms_coeff)
            @inbounds terms_coeff[i] *= -1.0
        end
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    else
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
    end
end

function Relax_Nonconvex_Quadratic!(trg, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                    lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64})
    terms_coeff = Float64[]
    terms_index = Int[]

    VarNum = length(n.LowerVar)
    quadratic_coefficient = zeros(Float64,VarNum)
    quadratic_constant = func.constant
    for term in func.quadratic_terms
        a = term.coefficient
        NodeIndx1 = src.VItoSto[term.variable_index_1.value]
        NodeIndx2 = src.VItoSto[term.variable_index_2.value]
        xL1 = n.LowerVar[NodeIndx1]
        xU1 = n.UpperVar[NodeIndx1]
        if NodeIndx1 == NodeIndx2               # quadratic terms
            if (a > 0.0)
                quadratic_coefficient[NodeIndx1] = a*(2.0*xL1+1.0)
                quadratic_constant -= a*xL1*xU1
            else
                quadratic_coefficient[NodeIndx1] = a*(xL1 + xU1)
                quadratic_constant -= a*xL1*xU1
            end
        else
            xL2 = n.LowerVar[NodeIndx2]
            xU2 = n.UpperVar[NodeIndx2]
            if (a > 0.0)
                check = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1] - xU1*xU2 + xL1*xL2 <= 0.0
                if check
                    quadratic_coefficient[NodeIndx1] += a*xL2
                    quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xL2
                else
                    quadratic_coefficient[NodeIndx1] += a*xU2
                    quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xU2
                end
            else
                check = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1] - xU1*xL2 + xL1*xU2 <= 0.0
                if check
                    quadratic_coefficient[NodeIndx1] += a*xL2
                    quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xL2
                else
                    quadratic_coefficient[NodeIndx1] += a*xU2
                    quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xU2
                end
            end
        end
    end

    # Adds affine terms
    for term in func.affine_terms
        quadratic_coefficient[src.VItoSto[term.variable_index.value]] += term.coefficient
    end

    temp = 0.0
    for (i,term) in enumerate(quadratic_coefficient)
        if (term != 0.0)
            push!(terms_coeff,term)
            push!(terms_index,src.StoToVI[i])
        end
    end

    varIndx = [MOI.VariableIndex(src.VItoSto[i]) for i in terms_index]

    if (lower == upper)
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,quadratic_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
        quadratic_constant *= -1.0
        for i in 1:length(terms_coeff)
            terms_coeff[i] *= -1.0
        end
        saf1_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf1 = MOI.ScalarAffineFunction{Float64}(saf1_term,quadratic_constant)
        MOI.add_constraint(trg,saf1,MOI.LessThan{Float64}(-lower))
    elseif (lower != -Inf)
        quadratic_constant *= -1.0
        for i in 1:length(terms_coeff)
            terms_coeff[i] *= -1.0
        end
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,quadratic_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    else
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,quadratic_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
    end
end

function RelaxQuadratic!(trg, src::Optimizer, n::NodeBB, r::RelaxationScheme)
    #count = 0
    x0 = (n.UpperVar - n.LowerVar)/2.0
    # Relax Convex Quadratic Terms
    for (func, set, ind) in src.QuadraticLEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(ind)
        if src.ConstraintConvexity[MOIindx]
            #count += 1; println("count: $count    LEQ Convex")
            Relax_Convex_Quadratic_Inner!(trg, src, func, -Inf, set.upper, n, x0)
        else
            #count += 1; println("count: $count    LEQ Nonconvex")
            Relax_Nonconvex_Quadratic!(trg, src, func, -Inf, set.upper, n, x0)
        end
    end
    for (func,set,ind) in src.QuadraticGEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(ind)
        if src.ConstraintConvexity[MOIindx]
            #count += 1; println("count: $count    GEQ Convex")
            Relax_Convex_Quadratic_Inner!(trg, src, func, set.lower, Inf, n, x0)
        else
            #count += 1; println("count: $count    GEQ Nonconvex")
            Relax_Nonconvex_Quadratic!(trg, src, func, set.lower, Inf, n, x0)
        end
    end
    for (func,set,ind) in src.QuadraticEQConstraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(ind)
        if src.ConstraintConvexity[MOIindx]
            #count += 1; println("count: $count    EQ Convex")
            Relax_Convex_Quadratic_Inner!(trg, src, func, set.value, set.value, n, x0)
        else
            #count += 1; println("count: $count    EQ Nonconvex")
            Relax_Nonconvex_Quadratic!(trg, src, func, set.value, set.value, n, x0)
        end
    end
end
