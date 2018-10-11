function UnivariateKernel(n::NodeBB,a::Float64,b::Float64,c::Float64,vi::Int)
        term1 = (c + b^2)/(4.0*a)
        term2 = term1/a
        if ((term1 > 0.0) && (a < 0.0)) # No solution, fathom node
            return false
        elseif (term2 => 0.0)
            xlo = n.LowerVar[vi]
            xhi = n.UpperVar[vi]
            chk1 = -sqrt(term2)-b/(2.0*a)
            chk2 = sqrt(term2)-b/(2.0*a)
            if (a > 0.0)
                (chk1 < xlo) && (n.LowerVar[vi] = max(xlo,chk2))
                (chk2 > xhi) && (n.UpperVar[vi] = min(xhi,chk1))
            else
                n.LowerVar[vi] = max(xlo,chk1)
                n.UpperVar[vi] = min(xhi,chk2)
        end
       (n.LowerVar[vi] <= n.UpperVar[vi])
end

function UnivariateQuadratic(m::Optimizer,n::NodeBB)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (a,b,c,vi) in m.UniQuadraticGEQConstraints
        feas = UnivariateKernel(n,a,b,c,vi)
        (~feas) && return feas
    end
    # fathom ax^2 + bx + c < u quadratics
    for (a,b,c,vi) in m.UniQuadraticLEQConstraints
         feas = UnivariateKernel(n,-a,-b,-c,vi)
         (~feas) && return feas
    end
    # fathom ax^2 + bx + c = v quadratics
    for (a,b,c,vi) in m.UniQuadraticEQConstraints
          feas = UnivariateKernel(n,a,b,c,vi)
          (~feas) && return feas
          feas = UnivariateKernel(n,-a,-b,-c,vi)
          (~feas) && return feas
     end
     return feas
end
#=
function BivariateKernel(n::NodeData,ax::Float64,ay::Float64,axy::Float64,
                         bx::Float64,by::Float64,vi::Int)
end

function BivariateQuadratic(m::Optimizer,n::NodeData)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (ax,ay,axy,bx,by,vi) in m.BiQuadraticGEQConstraints
        feas = BivariateKernel(n,ax,ay,axy,bx,by,vi)
        (~feas) && return feas
    end
    # fathom ax^2 + bx + c < u quadratics
    for (ax,ay,axy,bx,by,vi) in M.BiQuadraticLEQConstraints
         feas = BivariateKernel(n,ax,ay,axy,bx,by,vi)
         (~feas) && return feas
    end
    # fathom ax^2 + bx + c = v quadratics
    for (ax,ay,axy,bx,by,vi) in m.BiQuadraticEQConstraints
          feas = BivariateKernel(n,ax,ay,axy,bx,by,vi)
          (~feas) && return feas
          feas = BivariateKernel(n,ax,ay,axy,bx,by,vi)
          (~feas) && return feas
     end
     return feas
end


```
Checks to see if constraint is a bivariant quadratic term
```
function CheckBivariateQuad(f::MOI.ScalarQuadraticFunction{Float64})
    vIndx = Int64[]
    (length(f.quadratic_terms) > 3) && (return false)
    (length(f.affine_terms) > 2) && (return false)
    for i in f.affine_terms push!(vIndx,i.variable_index.value) end
    for i in f.quadratic_terms push!(vIndx,i.variable_index1.value, i.variable_index2.value) end
    (length(unique(vIndx)) == 2)
end

function SaveBivariate(func::MOI.ScalarQuadraticFunction{Float64},set::T) where {T<:MOI.AbstractScalarSet}
    vIndx = Int64[]
    for i in func.affine_terms push!(vIndx,i.variable_index.value) end
    for i in func.quadratic_terms push!(vIndx,i.variable_index1.value, i.variable_index2.value) end
    uIndx = unique(vIndx)
    vi1 = uIndx[1]
    vi2 = uIndx[2]
    ax =
    ay =
    axy =
    bx =
    by =
    c = func.constant
    ax,ay,axy,bx,by,c,vi
end

=#

```
Checks to see if constraint is a univariant quadratic term
```
function CheckUnivariateQuad(f::MOI.ScalarQuadraticFunction{Float64})
    (length(f.affine_terms) == 1) &&
    (length(f.quadratic_terms) == 1) &&
    (f.affine_terms[1].variable_index == f.quadratic_terms[1].variable_index_1 == f.quadratic_terms[1].variable_index_2)
end

GetValue(set::MOI.LessThan{Float64}) = set.upper
GetValue(set::MOI.GreaterThan{Float64}) = set.lower
GetValue(set::MOI.EqualTo{Float64}) = set.value

function GetUnivariateCoeff(func::MOI.ScalarQuadraticFunction{Float64},set::T) where {T<:MOI.AbstractScalarSet}
    a = func.quadratic_terms[1].coefficient
    b = (length(func.affine_terms) > 0) ?  func.affine_terms[1].coefficient : 0.0
    c = GetValue(set) - func.constant
    vi = func.quadratic_terms[1].variable_index.value
    a,b,c,vi
end

```
Classifies constraints as univariate or bivariate and adds them to storage vector
```
function ClassifyQuadratics!(m::Optimizer)
    # Check for Univariate and Bivariate Lesser Constraints
    for (func,set,indx) in QuadraticLEQConstraints
        if CheckUnivariateQuad(func)
            a,b,c,vi = GetUnivariateCoeff(func,set)
            push!(m.UniQuadraticLEQConstraints,(-a,-b,-c,VItoSto[vi]))
        #elseif CheckBivariateQuad(func)
        #    ax,ay,axy,bx,by,c,vi = GetBivariateCoeff(func,set)
        #    push!(BiQuadraticLEQConstraints,(ax,ay,axy,bx,by,c,VItoSto[vi]))
        end
    end

    # Check for Univariate and Bivariate Greater Constraints
    for (func,set,indx) in QuadraticGEQConstraints
        if CheckUnivariateQuad(func)
            a,b,c,vi = GetUnivariateCoeff(func,set)
            push!(UniQuadraticGEQConstraints,(a,b,c,VItoSto[vi]))
        #elseif CheckBivariateQuad(func)
        #    ax,ay,axy,bx,by,c,vi = GetBivariateCoeff(func,set)
        #    push!(BiQuadraticGEQConstraints,(ax,ay,axy,bx,by,c,VItoSto[vi]))
        end
    end

    # Check for Univariate and Bivariate Equality Constraints
    for (func,set,indx) in QuadraticEQConstraints
        if CheckUnivariateQuad(func)
            a,b,c,vi = GetUnivariateCoeff(func,set)
            push!(UniQuadraticEQConstraints,(a,b,c,VItoSto[vi]))
        #elseif CheckBivariateQuad(func)
        #    ax,ay,axy,bx,by,c,vi = GetBivariateCoeff(func,set)
        #    push!(BiQuadraticEQConstraints,(ax,ay,axy,bx,by,c,VItoSto[vi]))
        end
    end
end
