#' @name sample_lasso
#' @title Sample iid Lasso solutions in a non-uniqueness regime
#' @description
#' This function samples iid points from a set with non-unique elements,
#' \eqn{\Theta = \{ \theta\in\mathbb{R}^d: \|\theta\|_1 = \|\tilde\theta\|_1, \ X\tilde\theta = X\theta  \}}
#' whose strictly feasible representation is given by a object of class \code{sf_rep_lasso}, 
#' where \eqn{\tilde\theta} is a particular solution.
#' 
#' @param obj An S3 object of class \code{sf_rep_lasso}.
#' @param npoints The number of points that the function is going to sample.
#' @param random_walk Optional. A list that declares the random walk and some
#' related parameters. See the argument \code{random_walk} in \link[volesti]{sample_points} 
#' for details.
#' @param distribution Optional. A list that declares the target density and some 
#' related parameters. See the argument \code{distribution} in \link[volesti]{sample_points}
#'for details.
#' @return A matrix with \code{npoints} rows and \code{ncol(X)} columns. Each row 
#' is a sample point from the polytope \eqn{\Theta}.
#' 
#' @examples
#' # generate a sparse data set with non-unique Lasso solns
#' set.seed(1234)
#' n = 1000
#' d = 20
#' s = 5
#' rho = 0.01
#' Sigma = matrix(0, nrow = d, ncol = d)
#' for(i in 1:d) {
#'   for(j in 1:d) {
#'     Sigma[i,j] = rho^abs(i-j)
#'   }
#' }
#' X = mvtnorm::rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
#' X[,2] = -X[,1] 
#' X[,3] = X[,1]
#' beta = c(1, 0, 0, rep(1, s-3), rep(0, d-s))
#' epsilon = rnorm(n, mean = 0, sd = 0.1)
#' y = X %*% beta + epsilon

#' # find one particular solution 
#' lambda = 0.01
#' model = glmnet::glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, 
#'                        intercept = FALSE, standardize = FALSE)
#' beta_fit = as.numeric( model$beta )
#' 
#' # generate an S3 object of class sf_poly_lasso
#' equidx = equi_index_lasso(X, y, lambda, beta_fit, tol = 1e-3)
#' fit = sf_rep_lasso(X, beta_fit, equidx = equidx)
#' 
#' # sample 100 uniformly distributed points 
#' # from the predicted theta set
#' sample_lasso(fit, npoints = 100)
#' 
#' @export
sample_lasso = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  UseMethod("sample_lasso")
}

sample_lasso.default = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  stop("The input obj must be an S3 object of class sf_rep_lasso!")
}

sample_lasso.sf_rep_lasso = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  # do sampling
  P = volesti::Hpolytope(A = obj$Gamma, b = obj$gamma)
  samples = volesti::sample_points(P, n = npoints, random_walk = random_walk, distribution = distribution) 
  
  # convert back
  exe_convert = lasso_convert(samples, obj$Tm, obj$nv, npoints)
  exe_transtheta = lasso_trans_theta(exe_convert, obj$I, length(obj$equidx)) 
  result = matrix(1, nrow = npoints, ncol = 1) %*% t(as.matrix(obj$beta))
  result[, obj$equidx] = exe_transtheta
  return(result)
}
