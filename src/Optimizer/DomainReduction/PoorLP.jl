#Runs the Poor Man's LP to refine variables present nonlinear using linear constraints
function PoorManLP(m::Optimizer,n::NodeBB)

    # Runs Poor Man LP on constraints of form ax >= b
    for (func,constr,ind) in m.LinearGEQConstraints
        TempValue = (constr.lower - func.constant)
        for term in func.terms
            ti = m.VItoSto[term.variable_index.value]
            TempValue += -max(term.coefficient*n.UpperVar[ti], term.coefficient*n.LowerVar[ti])
        end
        for term in func.terms
            vi = m.VItoSto[term.variable_index.value]
            TermValue = -max(term.coefficient*n.UpperVar[vi], term.coefficient*n.LowerVar[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient > 0.0 )
                if (n.LowerVar[vi] < CutValue)
                    (CutValue > n.UpperVar[vi]) && (return false)
                    n.LowerVar[vi] = CutValue
                end
            else
                if (n.UpperVar[vi] > CutValue)
                    (CutValue < n.LowerVar[vi]) && (return false)
                    n.UpperVar[vi] = CutValue
                end
            end
        end
    end
        # Runs Poor Man LP on constraints of form ax <= b
    for (func,constr,ind) in m.LinearLEQConstraints
        TempValue = (constr.upper - func.constant)
        for term in func.terms
            ti = m.VItoSto[term.variable_index.value]
            TempValue += -min(term.coefficient*n.UpperVar[ti], term.coefficient*n.LowerVar[ti])
        end
        for term in func.terms
            vi = m.VItoSto[term.variable_index.value]
            TermValue = -min(term.coefficient*n.UpperVar[vi], term.coefficient*n.LowerVar[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient < 0.0 )
                if (n.LowerVar[vi] < CutValue)
                    (CutValue > n.UpperVar[vi]) && (return false)
                    n.LowerVar[vi] = CutValue
                end
            else
                if (n.UpperVar[vi] > CutValue)
                    (CutValue < n.LowerVar[vi]) && (return false)
                    n.UpperVar[vi] = CutValue
                end
            end
        end
    end

    for (func,constr,ind) in m.LinearEQConstraints
        TempValue = (constr.value - func.constant)
        for term in func.terms
            ti = m.VItoSto[term.variable_index.value]
            TempValue += -max(term.coefficient*n.UpperVar[ti], term.coefficient*n.LowerVar[ti])
        end
        for term in func.terms
            vi = m.VItoSto[term.variable_index.value]
            TermValue = -max(term.coefficient*n.UpperVar[vi], term.coefficient*n.LowerVar[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient > 0.0 )
                if (n.LowerVar[vi] < CutValue)
                    (CutValue > n.UpperVar[vi]) && (return false)
                    n.LowerVar[vi] = CutValue
                end
            else
                if (n.UpperVar[vi] > CutValue)
                    (CutValue < n.LowerVar[vi]) && (return false)
                    n.UpperVar[vi] = CutValue
                end
            end
        end
        TempValue = (constr.upper - func.constant)
        for term in func.terms
            ti = m.VItoSto[term.variable_index.value]
            TempValue += -min(term.coefficient*n.UpperVar[ti], term.coefficient*n.LowerVar[ti])
        end
        for term in func.terms
            vi = m.VItoSto[term.variable_index.value]
            TermValue = -min(term.coefficient*n.UpperVar[vi], term.coefficient*n.LowerVar[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient < 0.0 )
                if (n.LowerVar[vi] < CutValue)
                    (CutValue > n.UpperVar[vi]) && (return false)
                    n.LowerVar[vi] = CutValue
                end
            else
                if (n.UpperVar[vi] > CutValue)
                    (CutValue < n.LowerVar[vi]) && (return false)
                    n.UpperVar[vi] = CutValue
                end
            end
        end
    end

    return true
end
