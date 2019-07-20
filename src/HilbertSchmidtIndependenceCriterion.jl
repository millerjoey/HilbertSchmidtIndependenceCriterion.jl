module HilbertSchmidtIndependenceCriterion

  using StatsFuns, Distances

  # load Base modules
  using Statistics, LinearAlgebra
  import StatsBase: sample
	# includes
  include("common.jl")
  include("gammaHSIC.jl")

	# function exports
  export gammaHSIC, estimateKernelSize, rbfDotProduct

end # module
