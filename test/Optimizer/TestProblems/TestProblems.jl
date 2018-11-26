# Run LP tests
for i in 1:4
    include("LP/Prob$i.jl")
end

# Run QP tests
for i in 1:1
    include("QP/Prob$i.jl")
end

#=
# Run NLP tests
for i in 1:8
    include("NLP/Prob$i.jl")
end
=#
