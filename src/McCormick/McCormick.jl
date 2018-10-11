@reexport module McCormick

using IntervalArithmetic, StaticArrays, CommonSubexpressions, DiffRules, BenchmarkTools

import Base: +, -, *, /, convert, in, isempty, one, zero, real, eps, max, min,
             abs, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p, acosh,
             sqrt, sin, cos, tan, min, max, sec, csc, cot, ^, step, sign

import IntervalArithmetic: dist, mid, pow

# Export forward operators
export MC, cc, cv, Intv, lo, hi,  cc_grad, cv_grad, cnst, +, -, *, /, convert,
       one, zero, dist, real, eps, mid, exp, exp2, exp10, expm1, log, log2,
       log10, log1p, acosh, sqrt, sin, cos, tan, min, max, sec, csc, cot, ^,
       abs, step, sign, pow, in, isempty
       #acos, asin, atan, sinh, cosh, tanh, asinh, atanh, inv, sqr

# Export inplace operators
export plus!, mult!, min!, max!, minus!, div!, exp!, exp2!, exp10!, expm1!,
       log!, log2!, log10!, log1p!, sin!, cos!, tan!, asin!, acos!, atan!,
       sinh!, cosh!, tanh!, asinh!, acosh!, atanh!, abs!, sqr!, sqrt!, pow!

export seedg, IntervalType

# Export reverse operators
#=
export plus_rev!, mult_rev!, min_rev!, max_rev!, minus_rev!, div_rev!, exp_rev!,
       exp2_rev!, exp10_rev!, exp1m_rev!, log_rev!, log2_rev!, log10_rev!,
       log1p_rev!, sin_rev!, cos_rev!, tan_rev!, asin_rev!, acos_rev!, atan_rev!,
       sinh_rev!, cosh_rev!, tanh_rev!, asinh_rev!, acosh_rev!, atanh_rev!,
       abs_rev!, sqr_rev!, sqrt_rev!, pow_rev!
=#

# Export utility operators
#=
export grad, zgrad, ∩, mid3, MC_param, mid_grad, seed_g, line_seg, dline_seg,
       outer_rnd, cut, set_valid_check, set_subgrad_refine, set_multivar_refine,
       set_outer_rnd, tighten_subgrad, set_iterations, set_tolerance,
       set_diff_relax, default_options, value, mincv, maxcc, promote_rule
=#

include("ConvexityRules/ConvexityRules.jl")

include("Utilities/Constants.jl")
include("Utilities/InnerUtilities.jl")
include("Utilities/FastIntervals.jl")

include("OperatorLibrary/Type.jl")

include("Utilities/APIUtilities.jl")
include("Utilities/RootFinding.jl")

include("OperatorLibrary/ForwardOperators/Forward.jl")
#include("OperatorLibrary/ReverseOperators/Reverse.jl")
#include("OperatorLibrary/Inplace.jl")

end
