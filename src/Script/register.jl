function register(m::MOI.AbstractOptimizer, s::Symbol, dimension::Integer, f::Function)
    initNLP(m)
    if dimension == 1
        fprime = x -> ForwardDiff.derivative(f, x)
        fprimeprime = x -> ForwardDiff.derivative(fprime, x)
        Derivatives.register_univariate_operator!(m.nlp_data.user_operators, s, f, fprime, fprimeprime)
    else
        m.nlp_data.largest_user_input_dimension = max(m.nlp_data.largest_user_input_dimension,dimension)
        Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, s, UserAutoDiffEvaluator(dimension, f))
    end
end
