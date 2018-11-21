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
             gtemp = (x...) -> g([x[i] for i in 1:length(x)],pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp; autodiff=true)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(LBP_vars...) + $eps <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(LBP_vars...) + $eps))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.add_NL_constraint(jmodel, :(($gtemp_sym)(LBP_vars...) + $eps <= $tempU))
             end
         end
    return LBP_vars,jmodel
end

function loadscript_bnds_ubd!(m, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, f, g, pIndx, eps::Float64)

        # create JuMP model
        jmodel = Model(with_optimizer(m))

         # Add variables and associated bounds
         UBP_vars= @variable(jmodel, xL[i] <= UBP_vars[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f; autodiff=true)
         JuMP.set_NL_objective(jmodel, MOI.MinSense, Expr(:call, :f, UBP_vars...))

         println("finished objective")
         # Add nonlinear constraint
         for i in 1:ng
             gtemp = (x...) -> g([x[i] for i in 1:length(x)],pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp; autodiff=true)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 out = :($tempL <= ($gtemp_sym)(UBP_vars) + $eps <= $tempU)
                 println("out: $out")
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(UBP_vars) + $eps <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 out = :($tempL <= ($gtemp_sym)(UBP_vars) + $eps)
                 println("out: $out")
                 JuMP.add_NL_constraint(jmodel, :($tempL <= ($gtemp_sym)(UBP_vars) + $eps))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                out = :(($gtemp_sym)(UBP_vars) + $eps <= $tempU)
                println("out: $out")
                 JuMP.add_NL_constraint(jmodel, :(($gtemp_sym)(UBP_vars) + $eps <= $tempU))
             end
         end
    return UBP_vars,jmodel
end

function loadscript_llp!(m, np::Int, pL::Vector{Float64}, pU::Vector{Float64},
                         gin, xIndx, prob::Int)

         # create JuMP model
         jmodel = Model(with_optimizer(m))

         # Register nonlinear objective
         g(p...) = gin(xIndx,[p[i] for i in 1:length(p)])
         JuMP.register(jmodel, :g, np, g; autodiff=true)

         # Add variables and associated bounds
         if prob == 1
             var = @variable(jmodel, pL[i] <= LLP1_vars[i in 1:np] <= pU[i])
             JuMP.set_NL_objective(jmodel, MOI.MaxSense, Expr(:call, :g, LLP1_vars...))
             return LLP1_vars,jmodel
         else
             var = @variable(jmodel, pL[i] <= LLP2_vars[i in 1:np] <= pU[i])
             JuMP.set_NL_objective(jmodel, MOI.MaxSense, Expr(:call, :g, LLP2_vars...))
             return LLP2_vars,jmodel
         end
end
