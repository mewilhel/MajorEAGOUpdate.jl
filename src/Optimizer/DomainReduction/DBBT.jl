function Variable_DBBT!(x::NodeBB,
                        mult_lo::Vector{Float64},
                        mult_hi::Vector{Float64},
                        LBD::Float64,
                        UBD::Float64)
  nx::Int = length(x.LowerVar)
  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD <= UBD)
    for i = 1:nx
      if (mult_lo[i] > 0.0)
        lower_cut = x.UpperVar[i] - (UBD - LBD)/mult_lo[i]
        (lower_cut > x.LowerVar[i]) && (x.LowerVar[i] = lower_cut)
      elseif (mult_hi[i] > 0.0)
        upper_cut = x.LowerVar[i] + (UBD - LBD)/mult_lo[i]
        (upper_cut < x.UpperVar[i]) && (x.UpperVar[i] = upper_cut)
      end
    end
  end
  println("ran Variable_DBBT")
end
