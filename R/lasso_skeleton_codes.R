## equi-correlation ----
#' @name equi_index_lasso
#' @title The equi-correlation set of a Lasso problem
#' @description
#' This function finds the equi-correlation set of a Lasso problem. Given the design matrix 
#' \eqn{X\in\mathbb{R}^{n\times d}}, the response variable \eqn{y\in\mathbb{R}^n}, 
#' the tuning parameter \eqn{\lambda}, and the fitted coefficient \eqn{\hat\beta}, 
#' the equi-correlation set is defined as 
#' 
#' \eqn{\mathcal{E}= \{ i\in\{1,\dots, d\}: |X^\top_i(y - X\hat\beta)| = n\lambda \} }.
#' 
#' @details
#' The equi-correlation set includes indices for all relevant features. The detection of
#' \eqn{\mathcal{E}} depends on the accuracy of \eqn{\hat\beta} as the minimum point of
#' the Lasso loss function. The argument \code{tol} must be chosen to transcend the gap
#' between \eqn{n\lambda} and \eqn{|X^\top_i(y - X\hat\beta)|}. If the tolerance \code{tol}
#' is chosen too large, some irrelevant features may be included and thus the subsequent facial
#' reduction procedure would be decelerated; on the other hand, if the tolerance is chosen too 
#' small, some relevant features may be ignored which leads to a wrong strictly feasible representation. 
#' In conclusion, a too tight tolerance is more detrimental than a too loose one. An accurate
#' \eqn{\hat\beta} (\code{beta}) and a slightly loose \code{tol} are always helpful.
#' 
#' @param X the design matrix.
#' @param y a vector denoting the response variable.
#' @param lambda a numeric denoting the tuning parameter.
#' @param beta a vector denoting a particular Lasso solution.
#' @param tol a tolerance numeric greater than 0, which must cover the difference
#'  between \eqn{\lambda} and \eqn{|X^\top_i(y - X\hat\beta)|}.
#' 
#' @return An integer vector containing all elements of the equi-correlation set.
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
#' beta = as.numeric( model$beta )
#' 
#' # find the equi-correlation set
#' equi_index_lasso(X, y, lambda, beta, tol = 1e-5)
#' 
#' @export

equi_index_lasso = function(X, y, lambda, beta, tol = 1e-2) {
  # check
  if(!is.matrix(X)) {
    stop("X must be a matrix!")
  }
  if(!is.numeric(y)) {
    stop("y must be a numeric vector!")
  }
  if(!is.numeric(lambda) | length(lambda) != 1) {
    stop("lambda must be a length-one numeric vector!")
  }
  if(!is.numeric(beta)) {
    stop("beta must be a numeric vector!")
  }
  if(nrow(X) != length(y)) {
    stop("X and y are not comparable in size!")
  }
  if(ncol(X) != length(beta)) {
    stop("X and beta are not comparable in size!")
  }
  # begin
  n = nrow(X)
  z = abs(t(X) %*% (y - X %*% beta) / n)
  equal = sapply(z, function(j) {
    all.equal(lambda, j, tolerance = tol)
  }) 
  ind = seq_along(z)[equal == "TRUE"]
  return(ind)
}



gen_inv_index = function(X, beta, tol = 1e-2){
  b = X %*% beta
  z = pracma::pinv(X) %*% b
  equal = sapply(z, function(j) {
    all.equal(0, j, tolerance = tol)
  }) 
  ind = seq_along(z)[equal != "TRUE"]
  return(ind)
}



## right signs ----
check_sign = function(X, beta, n, d) {
  Pi = rbind(
    cbind(X, -X),
    cbind(t( rep(1, d) ) , t( rep(1, d) )),
    cbind(-X, X),
    cbind(t( rep(-1, d) ), t( rep(-1, d) )),
    diag(rep(-1, 2*d))
  )
  pi = c(
    X %*% beta,
    sum(abs(beta)),
    -X %*% beta,
    -sum(abs(beta)),
    rep(0, 2*d)
  )
  Omega = cbind( Pi, diag(2*(n+d+1)) )
  Omega = Omega[, -( (2*d+1):(2*d+2*n+2) )]
  I = lasso_fr(Omega, matrix(pi)) + 1 - 2 * d
  J = I + 2 * d
  K = I + 2 * n + 2
  Omega = Omega[-K , -union(I,J)]
  omega = pi[-K]
  return(list(
    I = I,
    Omega = Omega,
    omega = omega
  ))
}



report_sign = function(origd, equidx, I) {
  s = length(equidx)
  equi_sign = rep(1, 2 * s)
  equi_sign[I] = 0
  equi_sign = equi_sign[1:s] - equi_sign[(s+1):(2*s)]
  all_sign = rep(0, origd)
  all_sign[equidx] = equi_sign
  return(all_sign)
}



## intrinsic dim ----
lasso_int_dim = function(Omega, omega) {
  fullrank = ind_rows(Omega)
  return(list(
    N = fullrank$N,
    v = omega[fullrank$idx],
    indim = ncol(fullrank$N) - nrow(fullrank$N)
  ))
}



## minimal representation ----
lasso_polysf = function(N, v, I, d) {
  r = nrow(N)
  m = 2*d - length(I)
  E = solve(add_to_full(N))
  gamma = as.numeric( E[(m + 1):(2 * m), 1:r] %*% v )
  Gamma = as.matrix( -E[(m + 1):(2 * m), (r + 1):(2 * m)] )
  Tm = E[1:m, ]
  
  return(list(
    Gamma = Gamma,
    gamma = gamma,
    Tm = Tm
  ))
}



## convert ----
lasso_convert = function(samples, Tm, v, npoints) {
  nu = t( Tm %*% rbind( v %*% matrix(rep(1, npoints), nrow = 1), samples) )
  return(nu)
}



## trans_theta ----
lasso_trans_theta = function(nu, I, d) {
  Ic = setdiff(seq_len(2*d), I)
  Beta = matrix(0, nrow = nrow(nu), ncol = 2 * d)
  Beta[,Ic] = nu
  # Theta = sapply(1:d, function(j) {
  #   Beta[,j] - Beta[,j+d]
  # })
  Beta = Beta[,1:d] - Beta[,(d+1):(2*d)]
  return(Beta)
}
