#' @name sf_rep
#' @title Strictly feasible representation of a polytope
#' @description
#' This function finds a strictly feasible representation of the given polytope 
#' \eqn{\mathcal{S} = \{x\in\mathbb{R}^n:Ax\leq b\}} within a subspace of its intrinsic dimension.
#' 
#' @details
#' Given a polytope \eqn{\mathcal{S} = \{x\in\mathbb{R}^n:Ax\leq b\}}, namely a 
#' bounded intersection of finite number of half-spaces where each half-space is 
#' represented by a constraint \eqn{\{x:a_i^\top x\leq b_i\}}. A polytope 
#' \eqn{\mathcal{S}} with implicit equality constraints is of measure 0 in \eqn{\mathbb{R}^n}
#' and has no strictly feasible points (or interior points) inside \eqn{\mathcal{S}}. 
#' This degeneration may fail standard sampling algorithms over a polytope, including the
#' Ball Walk, the Hit-and-Run, and the Dikin Walk.
#' 
#' The function \code{sf_rep} finds a strictly feasible representation of any polytope \eqn{\mathcal{S}} 
#' with at least one explicit inequality constraint. Specifically, this function 
#' identifies the implicit equality constraints of \eqn{\mathcal{S}}, 
#' presents its intrinsic dimension \eqn{{\rm dim}\mathcal{S}}, 
#' and provides an equivalent strictly feasible polytope 
#' \eqn{\mathcal{S}^* = \{y\in\mathbb{R}^{{\rm dim}\mathcal{S}}:\Gamma y\leq \gamma\}}.
#' In addition, a transformation matrix pair \eqn{(T, d)} is given to convert \eqn{\mathcal{S}^*}
#' back to  \eqn{\mathcal{S}}: For any \eqn{y\in\mathcal{S}^*}, we have 
#' \eqn{T \begin{bmatrix} d \\ y \end{bmatrix} \in \mathcal{S}}.
#' 
#' @param A The matrix \eqn{A} of the polytope \eqn{\mathcal{S}}.
#' @param b The matrix \eqn{b} of the polytope \eqn{\mathcal{S}}.
#' @return The function \code{sf_rep} returns a S3 object of class sf_rep
#'  containing the following components:
#'   \item{A}{The argument \code{A}.}
#'   \item{b}{The argument \code{b}.}
#'   \item{I}{An integer vector containing the indices of implicit constraints.}
#'   \item{indim}{An integer denoting the intrinsic dimension of \eqn{\mathcal{S}}.}
#'   \item{Gamma}{The matrix \eqn{\Gamma} of the polytope \eqn{\mathcal{S}^*}. }
#'   \item{gamma}{The vector \eqn{\gamma} of the polytope \eqn{\mathcal{S}^*}. }
#'   \item{Tm}{The matrix \eqn{T} in the transformation pair \eqn{(T,d)}.}
#'   \item{d}{The vector \eqn{d} in the transformation pair \eqn{(T,d)}.}

#' @examples
#' # generate A and b of the polytope
#' A = matrix(c(1, 1, -1, -1, -1, 0, 0, -1), byrow = TRUE, nrow = 4)
#' b = c(1, -1, 0, 0)
#' 
#' # find the strictly feasible representation
#' sf_rep(A, b)
#' @export

sf_rep = function(A, b) {
  # preparation
  m = nrow(A)
  n = ncol(A)
  # implicit equality
  I = im_equ(A, matrix(b)) + 1
  M = cbind(A, diag(m))[,-I]
  I = I - n
  if( length(I) == 0 ) {
    stop("The polytope has no implicit equality constraints.")
  }
  if( length(I) == m) {
    stop("All constraints are implicit equalities.")
  }
  
  # intrinsic dimension
  fullrank = ind_rows(M)
  N = fullrank$N
  d = b[fullrank$idx]
  r = nrow(N)
  # polysf
  E = solve( add_to_full(N) )
  Gamma = matrix( -E[(n+1):ncol(N), (r+1):ncol(N)], nrow = ncol(N) - n )
  gamma = matrix( E[(n+1):ncol(N), 1:r], nrow = ncol(N) - n ) %*% d
  Tm = E[1:n,]
  
  ls = list(
    A = A,
    b = b,
    I = I,
    indim = ncol(N) - nrow(N),
    Gamma = Gamma,
    gamma = as.vector( gamma ),
    Tm = Tm,
    d = d
  )
  class(ls) = "sf_rep"
  return(ls)
}


print.sf_rep = function(x, ...) {
  cat("Indices for implicit equality constraints (rows): ","\n")
  cat(x$I)
  cat("\n\n")
  
  cat("Intrinsic dimension: ","\n")
  cat(x$indim)
  cat("\n\n")
}
