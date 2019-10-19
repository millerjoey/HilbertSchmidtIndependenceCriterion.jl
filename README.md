# Hilbert-Schmidt Independence Criterion (HSIC)

This package provides basic implementations of the Hilbert-Schmidt Independence Criterion (HSIC) for Julia 1.0.

## What is implemented
The package currently contains the following implementations:

- Gamma HSIC (HSIC with Gamma approximation) [1]

## Example

The gamma HSIC can be run using:

```julia
X = randn(100) # rows are samples
Y = randn(100) * 0.2 # rows are samples
p = gammaHSIC(X, Y)
threshold = 0.05
independent = p < threshold
```

## Reference
[1] Gretton, Arthur, et al. "A kernel statistical test of independence." Advances in Neural Information Processing Systems. 2007.
