# Following Case Fails: sin(x[1]) + y * 3.0*2.0*x[2]
# Next Case: Assignment in Script
# Next Case: For loop
# Next Case: While loop
# Next Case: If-else
# Next Case: User-defined...

module Tracer

    import Base: abs, sin, cos, tan, sec, csc, cot, asin, acos, atan, asec, acsc,
              acot, sinh, cosh, tanh, asinh, acosh, atanh, sech, asech, csch,
              acsch, coth, acoth, sqrt, log, log2, log10, log1p, exp, exp2,
              exp10, expm1, +, -, step, sign, real, inv, *,
              +, -, / , min, max, >, <, ==, =>, <= ,^, getindex, afoldl

    using Cassette, JuMP, SparseArrays

    include("src/Types.jl")
    include("src/Utilities.jl")

    Cassette.@context TraceCtx

    include("src/ExecuteOps.jl")
    include("src/ExecuteUtils.jl")

    # getindex(x,1) to getindex(x,n) correspond to 1..n variable
    function trace_script(f,n)
        tape = Tape(n)
        x = SetTraceSto([SetTrace(i) for i=1:n])
        Cassette.overdub(TraceCtx(metadata = tape), f, x)
        return tape
    end

    include("src/ToJuMP.jl")

end # module

using JuMP, Ipopt

function f(x)
    #z = abs(x[1]) + 3.0*x[2]
    y = 1.3
    #for i in 1:3
    #    y += i*z
    #end
    #cos(x[1]) + y
    return sin(3.0*x[1]) + y*3.0*x[2]
end

tape = Tracer.trace_script(f,2)
Tracer.child_to_parent!(tape)

m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m,  0 <= x[1:2] <= 1)
@variable(m, y)
@NLobjective(m, Min, sin(3.0*x[1]) + 3.9*x[2])
JuMP.optimize!(m)

moi_backend = JuMP.backend(m)
evaluator = moi_backend.optimizer.model.optimizer.nlp_data.evaluator
objstorage =  evaluator.objective
#+ 3.0*x[2]
