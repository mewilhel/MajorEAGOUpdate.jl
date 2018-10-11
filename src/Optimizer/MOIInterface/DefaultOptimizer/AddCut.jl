function CheckValidSafety!()
end

function IsSafeLinearCut(x::Optimizer,a::Vector{Float64},b::Float64)
    (abs(b) > x.CutSafetyConstant) && (return false)
    lena = length(a)
    for i in 1:lena
        if a[i] != 0.0
            ~(x.CutSafetyLower <= abs(a[i]) <= x.CutSafetyUpper) && (return false)
            for j in 1:lena
                if a[j] != 0.0
                    ~(x.CutSafetyLower <= abs(a[i]/a[j]) <= x.CutSafetyUpper) && (return false)
                end
            end
        end
    end
    return true
end


function DefaultAddCut(x::Optimizer)
end
