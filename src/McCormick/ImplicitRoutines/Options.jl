"""
    mc_opts

Storage type used for parameter options in implicit bounding routine.
* `lambda::Float64`: Affine weighting parameter
* `kmax::Int`: Number of iterations run
* `LAlg::Symbol`:
* `CTyp::Symbol`:
* `np::Int`:
* `nx::Int`:
* `aff_correct_eps::Float64`: Affine correction tolerance
"""
mutable struct mc_opts
  lambda::Float64
  kmax::Int
  LAlg::Symbol   # Type of symbol
  CTyp::Symbol   #
  np::Int
  nx::Int
  aff_correct_eps::Float64
end

"""
    mc_opts()

Initialization function for currently sets the weight `.lambda = 0.5`, the
contractor style `.style = KrawczykCW`, the number of iterations to `.kmax = 2`,
and the affine correction tolerance as `.aff_correct_eps = 1E-12`. Other rounding
options disabled.
"""
mc_opts() = mc_opts(0.5,1,:Dense,:Newton,0,0,1E-10)

"""
    set_default!(x::mc_opts)

Restores `x` to default values. Sets the weight `.lambda = 0.5`, the
contractor style `.style = KrawczykCW`, the number of iterations to `.kmax = 2`,
and the affine correction tolerance as `.aff_correct_eps = 1E-12`. Other rounding
options disabled.
"""
function set_default!(x::mc_opts)
  x.lambda = 0.5
  x.kmax = 2
  LAlg = :Dense
  CTyp = :Krawczyk
#x.aff_correct_eps = 1E-12
end
