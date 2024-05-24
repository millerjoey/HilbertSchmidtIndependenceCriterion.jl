"""
    estimateKernelSize(X::AbstractArray; sampleSize = 100)

    Estimate kernel size by taking median distance between points.
"""
function estimateKernelSize(X::AbstractArray; sampleSize = 100)

	M = length(X)

	# set kernel size to median distance between points
	if M > sampleSize
		Xmed = reshape(sample(X, sampleSize, replace=false), (sampleSize, 1))
		S = sampleSize
	else
		Xmed = reshape(X, (M, 1))
		S = M
	end

	dists = pairwise(SqEuclidean(), Xmed, Xmed, dims = 1)
    sig = sqrt(0.5 * median(dists))

    return sig
end

"""
    rbfDotProduct(X::Array, X::Array, kernelSize::Float64)

"""
function rbfDotProduct(X::Array, Y::Array, kernelSize)

	G = sum((X.*Y), dims = 2)

	Q = repeat(G, 1, length(Y))
	R = repeat(G', length(X), 1)

	H = Q + R - 2*X*Y'

	return exp.(-H/2/kernelSize^2)

end

eye(n) = Matrix(I, n, n)
