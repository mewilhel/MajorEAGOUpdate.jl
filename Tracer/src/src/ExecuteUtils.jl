# define primitives for associative terms
Cassette.execute(ctx::TraceCtx, ::typeof(*), x, y) = afoldl(x, y)
Cassette.execute(ctx::TraceCtx, ::typeof(afoldl), x, y, z) = afoldl(x, y, z)
Cassette.execute(ctx::TraceCtx, ::typeof(afoldl), a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...) = afoldl(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...)

# primitive for array access
Cassette.execute(ctx::TraceCtx, ::typeof(getindex), A::Array, i::Int) = getindex(A,i)
Cassette.execute(ctx::TraceCtx, ::typeof(getindex), A::SetTraceSto, i::Int) = getindex(A,i)

function Cassette.prehook(::TraceCtx, f, args...)
    println(f, args)
end
