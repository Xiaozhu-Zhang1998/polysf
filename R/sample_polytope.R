#' @name sample_polytope
#' @title Sample uniformly or normally distributed points from a polytope
#' @description
#' This function samples uniformly or normally distributed points
#' from a polytope \eqn{\mathcal{S} = \{x\in\mathbb{R}^n:Ax\leq b\}} whose 
#' strictly feasible representation is given by a object of class \code{sf_rep}.
#' @param obj An S3 object of class \code{sf_rep}.
#' @param npoints The number of points that the function is going to sample.
#' @param random_walk Optional. A list that declares the random walk and some
#' related parameters. See the argument \code{random_walk} in \link[volesti]{sample_points} 
#' for details.
#' @param distribution Optional. A list that declares the target density and some 
#' related parameters. See the argument \code{distribution} in \link[volesti]{sample_points}
#'for details.
#' @return A matrix with \code{npoints} rows and \code{ncol(A)} columns. Each row 
#' is a sample point from the polytope \eqn{\mathcal{S}}.
#' @examples 
#' # generate an S3 object of class sf_poly
#' A = matrix(c(1, 1, -1, -1, -1, 0, 0, -1), byrow = TRUE, nrow = 4)
#' b = c(1, -1, 0, 0)
#' fit = sf_rep(A, b)
#' 
#' # sample 100 uniformly distributed points from Ax<=b
#' sample_polytope(fit, npoints = 100)
#' 
#' @export
sample_polytope = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  UseMethod("sample_polytope")
}

sample_polytope.default = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  stop("The input obj must be an S3 object of class sf_rep!")
}

sample_polytope.sf_rep = function(obj, npoints, random_walk = NULL, distribution = NULL) {
  # do sampling
  P = volesti::Hpolytope(A = obj$Gamma, b = obj$gamma)
  samples = volesti::sample_points(P, n = npoints, random_walk = random_walk)
  nu = t( obj$Tm %*% rbind( matrix(obj$d) %*% matrix(rep(1, npoints), nrow = 1), samples) )
  return(nu)
}
