#' @name sf_rep_lasso
#' @title Strictly feasible representation of Lasso in a non-uniqueness regime
#' @description
#' Given the design matrix \eqn{X\in\mathbb{R}^{n\times d}} and a particular solution \eqn{\tilde\theta}, 
#' this function finds a strictly feasible representation of the set
#' \eqn{\Theta = \{ \theta\in\mathbb{R}^d: \|\theta\|_1 = \|\tilde\theta\|_1, \ X\tilde\theta = X\theta  \}}
#' which can be proven a polytope with implicit equality constraints.
#' 
#' @details
#' The function \code{sf_rep_lasso} can deal with two types of theta set:
#' \enumerate{
#'    \item The predicted theta set:
#' 
#'    \eqn{\hat\Theta = \{ \theta\in\mathbb{R}^d: \theta \in \arg\min_{\vartheta} \frac{1}{2n} \|X \vartheta - y\|_2^2 + \lambda \|\vartheta\|_1  \}};
#'    
#'    It can be shown that \eqn{\hat\Theta} is the same as \eqn{\Theta} when the particular solution 
#'    \eqn{\tilde\theta \in \hat\Theta}.
#'    
#'    \item The true theta set: given a particular true coefficient \eqn{\theta^*},
#' 
#'    \eqn{\Theta^* = \{ \theta\in\mathbb{R}^d: \|\theta\|_1 = \|\theta^*\|_1, \ X\theta = X\theta^*  \}}.
#' }
#' We use the generic notation \eqn{\Theta} with the particular solution \eqn{\tilde\theta} to denote both types.
#' Under certain conditions, \eqn{\Theta} may contain non-unique (and infinitely many) elements.
#' For both types the steps to obtaining their strictly feasible representation are exactly the same. 
#' However, there are some caveats in terms of determining the equi-correlation set:
#' \eqn{(X, \tilde\beta, y, \lambda)} are required for the predicted theta set,
#' while only \eqn{(X, \tilde\beta)} are needed for the true theta set. 
#' 
#' Essentially, the strictly feasible representation of \eqn{\Theta} is obtained through facial
#' reduction restricted to the equi-correlation set, but the procedure can be even simplified
#' due to the special structure of \eqn{(\Theta, \tilde\theta)}.
#' The result of facial reduction reveals the sign of the theta set which is given by
#' 
#' \eqn{{\rm sgn}(\Theta)_j = \begin{cases} 1,  & \theta_j > 0\ \exists \theta\in\Theta ,\\
#' -1, & \theta_j < 0\ \exists \theta\in\Theta ,\\
#' 0,  & \theta_j = 0\ \forall \theta\in\Theta .
#' \end{cases}} 
#' 
#' Note that all points of \eqn{\Theta} must be placed in same the orthants so \eqn{{\rm sgn}(\Theta)} 
#' is well-defined.
#' 
#' The strictly feasible representation of \eqn{\Theta} is given by
#' \eqn{\mathcal{K} = \{y\in\mathbb{R}^{{\rm dim}\mathcal{\Theta}}:\Gamma y\leq \gamma\}}
#' where \eqn{{\rm dim}\mathcal{\Theta}} is the intrinsic dimension of \eqn{\Theta}. Note that
#' there exists a bijection between \eqn{\Theta} and \eqn{\mathcal{K}} based on the transformation
#' pair \eqn{(T, \nu)}. Particularly, for any \eqn{y\in\mathcal{\Theta}}, we have 
#' \eqn{ T \begin{bmatrix} \nu \\ y \end{bmatrix} \in\Theta }.
#' 
#' 
#' 
#' @param X The design matrix.
#' @param beta A vector denoting a particular solution.
#' @param type A character, \code{"hat"} if the set \eqn{\Theta} 
#' is the predicted theta set; \code{"star"} if the set \eqn{\Theta} is the true theta set; 
#' can be left \code{NULL} if \code{equidx} is provided.
#' @param y A vector denoting the response variable. Need be provided only when 
#' \code{type = "hat"}.
#' @param lambda A numeric denoting the tuning parameter. 
#' Need be provided only when \code{type = "hat"}.
#' @param equidx Optional. An integer vector containing all elements of the equi-correlation set.
#' @param tol A tolerance numeric greater than 0. Need be provided only when \code{type = "hat"} 
#' or \code{type = "star"}. See \code{tol} in \link[polysf]{equi_index_lasso} for details.
#' 
#' @return The function \code{sf_rep_lasso} returns an S3 object of class \code{sf_rep_lasso}
#' containing the following components:
#'    \item{X}{The argument \code{X}.}
#'    \item{beta}{The argument \code{beta}.}
#'    \item{equidx}{The equi-correlation set either from input or found based on \code{type}, \code{y},
#'      \code{lambda}, and \code{tol}.}
#'    \item{I}{An integer vector for possible future sampling.}
#'    \item{signs}{An integer vector of length \eqn{d} denoting the signs of \eqn{\Theta}. For each feature,
#'      +1 means the feature is active and its coefficient is always non-negative; 0 means the feature is inactive;
#'      -1 means the feature is active and its coefficient is always non-positive.   }
#'    \item{indim}{An integer denoting the intrinsic dimension of \eqn{\Theta}. }
#'    \item{Gamma}{The matrix \eqn{\Gamma} of the polytope \eqn{\mathcal{K}}.}
#'    \item{gamma}{The matrix \eqn{\gamma} of the polytope \eqn{\mathcal{K}}.}
#'    \item{Tm}{The matrix \eqn{T} in the transofrmation pair \eqn{(T, \nu)}.}
#'    \item{nv}{The matrix \eqn{\nu} in the transofrmation pair \eqn{(T, \nu)}.}
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
#' # Ex1: find the SF Rep of the predicted theta set 
#' # with pre-calculated equi-correlation set
#' equidx = equi_index_lasso(X, y, lambda, beta_fit, tol = 1e-3)
#' sf_rep_lasso(X, beta_fit, equidx = equidx)
#' 
#' # without pre-calculated equi-correlation set
#' sf_rep_lasso(X, beta_fit, type = "hat", y = y, lambda = lambda, tol = 1e-3)
#' 
#' # Ex2: find the SF Rep of the true theta set 
#' sf_rep_lasso(X, beta, type = "star", tol = 1e-5)
#' 
#' 
sf_rep_lasso = function(X, beta, type = NULL, y = NULL, lambda = NULL, tol = 1e-2, equidx = NULL) {
  # save the beta
  beta_store = beta
  X_store = X
  
  # pre-check & find the equidx
  if(is.null(equidx)) {
    if(type == "hat") {
      if(is.null(y) | is.null(lambda)) {
        stop("y or lambda is missing!")
      }
      equidx = equi_index_lasso(X, y, lambda, beta, tol = tol)
    } else if(type == "star") {
      equidx = gen_inv_index(X, beta, tol = tol)
    } else {
      stop("The type is invalid!")
    }
  }
  
  # check whether all elements of equidx are 0
  if(length(equidx) == 0) {
    stop("The equi-correlation set is empty!")
  }
  
  # subset according to equi-correlation / gen_inv index
  X = as.matrix(X[, equidx])
  beta = beta[equidx]
  
  # go through the sign of each var
  n = nrow(X)
  d = ncol(X)
  exe_checksign = check_sign(X, beta, n, d)
  all_sign = report_sign(ncol(X_store), equidx, exe_checksign$I)
  
  # find intrinsic dimension & strict-feasible representation
  exe_intdim = lasso_int_dim(exe_checksign$Omega, exe_checksign$omega)
  if (exe_intdim$indim == 0) {
    stop("The LASSO solution is unique!")
  }
  exe_polysf = lasso_polysf(exe_intdim$N, exe_intdim$v, exe_checksign$I, d)
  
  # return object
  ls = list(
    X = X_store,
    beta = beta_store,
    equidx = equidx,
    I = exe_checksign$I,
    signs = all_sign,
    indim = exe_intdim$indim,
    Gamma = exe_polysf$Gamma,
    gamma = exe_polysf$gamma,
    Tm = exe_polysf$Tm,
    nv = exe_intdim$v
  )
  class(ls) = "sf_rep_lasso"
  return(ls)
}



print.sf_rep_lasso = function(x, ...) {
  cat("Equi-correlation set:", "\n")
  cat(x$equidx)
  cat("\n\n")
  
  cat("Sign of the Theta set:", "\n")
  cat(x$signs)
  cat("\n\n")
  
  cat("Intrinsic dimension:", "\n")
  cat(x$indim)
  cat("\n\n")
}
