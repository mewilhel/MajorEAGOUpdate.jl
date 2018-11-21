# Main decisions (replace nlp evaluator or recode the original EAGO/JuMP ones)

function loadscript_bnds_lbd!(m, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64},f, g, pIndx, LBP_vars)

        println("typeof(m): $(typeof(m))")

        # create JuMP model
        jmodel = Model(with_optimizer(m))

         # Add variables and associated bounds
         LBP_vars = @variable(jmodel, xL[i] <= LBP_vars[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f; autodiff=true)
         JuMP.set_NL_objective(jmodel, MOI.MinSense, Expr(:call, :f, LBP_vars...))

         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)(LBP_vars) + $eps <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)(LBP_vars) + $eps))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :(($gtemp_sym)(LBP_vars) + $eps <= $tempU))
             end
         end
    return LBP_vars,jmodel
end

function loadscript_bnds_ubd!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function, pIndx, eps::Float64)

        # create JuMP model
        jmodel = Model(m)

         # Add variables and associated bounds
         UBP_vars= @variable(jmodel, xL[i] <= UBP_vars[i in 1:nx] <= xU[i])

         # Add nonlinear objective (nothing defined with 3 arguments)
         JuMP.register(jmodel, :f, nx, f; autodiff=true)
         JuMP.set_NL_objective(jmodel, MOI.MinSense, Expr(:call, :f, UBP_vars...))

         println("finished objective")
         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(UBP_vars) + $eps <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(UBP_vars) + $eps))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.add_NL_constraint(jmodel, :(($gtemp_sym)(UBP_vars) + $eps <= $tempU))
             end
         end
    return var,jmodel
end

function loadscript_llp!(m::MOI.AbstractOptimizer, np::Int, pL::Vector{Float64}, pU::Vector{Float64},
                         g::Function, xIndx, sense::MOI.OptimizationSense, prob::Int)

         # create JuMP model
         jmodel = Model(m)

         # Register nonlinear objective
         JuMP.register(jmodel, :g, np, p -> g(xIndx,p))

         # Add variables and associated bounds
         if prob == 1
             var = @variable(jmodel, pL[i] <= LLP1_vars[i in 1:nx] <= pU[i])
             @NLobjective(jm, Min, g(LLP1_vars))
         else
             var = @variable(jmodel, pL[i] <= LLP2_vars[i in 1:nx] <= pU[i])
             @NLobjective(jm, Min, g(LLP2_vars))
         end

         return var,jmodel
end
