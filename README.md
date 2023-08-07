
# polysf

<!-- badges: start -->
<!-- badges: end -->

The package `polysf` is designed for providing a strictly feasible representation for a polytope with implicit equality constraints, and a sampler on the degenerated polytope. The same functionality for Lasso solutions in a non-uniqueness regime, as a special case of the degenerated polytope, is provided as well.

## Installation

``` r
devtools::install_github(repo = "Xiaozhu-Zhang1998/Polysf")
library(Polysf)
```

## Example

``` r
library(polysf)
```
### A degenerated polytope

``` r
# generate an S3 object of class sf_poly
A = matrix(c(1, 1, -1, -1, -1, 0, 0, -1), byrow = TRUE, nrow = 4)
b = c(1, -1, 0, 0)
fit = sf_rep(A, b)

# sample 100 uniformly distributed points from Ax<=b
sample_polytope(fit, npoints = 100)
```

### Lasso in a non-uniqueness regime

``` r
# generate a sparse data set with non-unique Lasso solns
set.seed(1234)
n = 1000
d = 20
s = 5
rho = 0.01
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = mvtnorm::rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
X[,2] = -X[,1]
X[,3] = X[,1]
beta = c(1, 0, 0, rep(1, s-3), rep(0, d-s))
epsilon = rnorm(n, mean = 0, sd = 0.1)
y = X %*% beta + epsilon
# find one particular solution 
lambda = 0.01
model = glmnet::glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda,
                       intercept = FALSE, standardize = FALSE)
beta_fit = as.numeric( model$beta )

# generate an S3 object of class sf_poly_lasso
equidx = equi_index_lasso(X, y, lambda, beta_fit, tol = 1e-3)
fit = sf_rep_lasso(X, beta_fit, equidx = equidx)

# sample 100 uniformly distributed points 
# from the predicted theta set
sample_lasso(fit, npoints = 100)
```
