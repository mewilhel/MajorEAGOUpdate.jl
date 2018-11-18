# Main decisions (replace nlp evaluator or recode the original EAGO/JuMP ones)

function loadscript_lbd!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function, pIndx, varnam)

        # create JuMP model
        jmodel = Model(m)

         # Add variables and associated bounds
         @variable(jmodel, xL[i] <= x[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f)
         @NLobjective(jm, Min, f(x))

         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam) <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam)))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :(($gtemp_sym)($varnam) <= $tempU))
             end

         end
end

function loadscript_llp1!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function, pIndx, varnam)

         # create JuMP model
         jmodel = Model(m)

         # Add variables and associated bounds
         @variable(jmodel, xL[i] <= x[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f)
         @NLobjective(jm, Min, f(x))

         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam) <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam)))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :(($gtemp_sym)($varnam) <= $tempU))
             end

         end
end

function loadscript_llp2!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function, pIndx, varnam)

         # create JuMP model
         jmodel = Model(m)

         # Add variables and associated bounds
         @variable(jmodel, xL[i] <= x[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f)
         @NLobjective(jm, Min, f(x))

         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam) <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam)))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :(($gtemp_sym)($varnam) <= $tempU))
             end

         end
end

function loadscript_ubd!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function, pIndx, varnam)

        # create JuMP model
        jmodel = Model(m)

         # Add variables and associated bounds
         @variable(jmodel, xL[i] <= x[i in 1:nx] <= xU[i])

         # Add nonlinear objective
         JuMP.register(jmodel, :f, nx, f)
         @NLobjective(jmodel, Min, f(x))

         # Add nonlinear constraint
         for i in 1:ng
             gtemp = x -> g(x,pIndx)            # ONLY FOR LOWER AND UPPER
             gtemp_sym = Symbol("g$i")
             JuMP.register(jmodel, gtemp_sym, nx, gtemp)
             if ((gL[i] != -Inf) && (gU[i] != Inf))
                 tempL = gL[i]; tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam) <= $tempU))
             elseif (gL[i] != -Inf)
                 tempL = gL[i];
                 JuMP.addNLconstraint(jmodel, :($tempL <= ($gtemp_sym)($varnam)))
             elseif (gU[i] != Inf)
                tempU = gU[i];
                 JuMP.addNLconstraint(jmodel, :(($gtemp_sym)($varnam) <= $tempU))
             end

         end
end

function loadscript_implicit!(m::MOI.AbstractOptimizer, nx::Int, ng::Int, xL::Vector{Float64}, xU::Vector{Float64},
                     gL::Vector{Float64}, gU::Vector{Float64}, sense::MOI.OptimizationSense,
                     f::Function, g::Function)
end
