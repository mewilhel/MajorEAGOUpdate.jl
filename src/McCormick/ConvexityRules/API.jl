const CV_TRAITS = Dict{Tuple{Union{Expr,Symbol},Symbol,Int},Tuple{Symbol,Symbol,Union{Symbol,Float64},Float64,Float64}}()

cvtrait(M::Union{Expr,Symbol}, f::Symbol, q::Symbol, args...) = CV_TRAITS[M,f,q,length(args)](args...)

function _get_quoted_symbol(ex::Expr)
    @assert ex.head == :quote
    @assert length(ex.args) == 1 && isa(ex.args[1], Symbol) "Function not a single symbol"
    ex.args[1]
end

function _get_quoted_symbol(ex::QuoteNode)
    @assert isa(ex.value, Symbol) "Function not a single symbol"
    ex.value
end

macro set_cvtrait(def,c,m,val,lb,ub)
    @assert isa(def, Expr)
    lhs = def
    @assert isa(lhs, Expr) && lhs.head == :call "LHS is not a dot call"
    qualified_f = lhs.args[1]
    @assert isa(qualified_f, Expr) && qualified_f.head == :(.) "Function is not qualified by module"
    M = qualified_f.args[1]
    f = _get_quoted_symbol(qualified_f.args[2])
    args = lhs.args[2:end]
    key = Expr(:tuple, Expr(:quote, M), Expr(:quote, f), length(args))
    rule = Expr(:tuple, c, m, val, lb, ub)
    return esc(quote
        $ConvexityRules.CV_TRAITS[$key] = $rule
        $key
    end)
end
