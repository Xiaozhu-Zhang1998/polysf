% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lasso_skeleton_codes.R
\name{equi_index_lasso}
\alias{equi_index_lasso}
\title{The equi-correlation set of a Lasso problem}
\usage{
equi_index_lasso(X, y, lambda, beta, tol = 0.01)
}
\arguments{
\item{X}{the design matrix.}

\item{y}{a vector denoting the response variable.}

\item{lambda}{a numeric denoting the tuning parameter.}

\item{beta}{a vector denoting a particular Lasso solution.}

\item{tol}{a tolerance numeric greater than 0, which must cover the difference
between \eqn{\lambda} and \eqn{|X^\top_i(y - X\hat\beta)|}.}
}
\value{
An integer vector containing all elements of the equi-correlation set.
}
\description{
This function finds the equi-correlation set of a Lasso problem. Given the design matrix 
\eqn{X\in\mathbb{R}^{n\times d}}, the response variable \eqn{y\in\mathbb{R}^n}, 
the tuning parameter \eqn{\lambda}, and the fitted coefficient \eqn{\hat\beta}, 
the equi-correlation set is defined as 

\eqn{\mathcal{E}= \{ i\in\{1,\dots, d\}: |X^\top_i(y - X\hat\beta)| = n\lambda \} }.
}
\details{
The equi-correlation set includes indices for all relevant features. The detection of
\eqn{\mathcal{E}} depends on the accuracy of \eqn{\hat\beta} as the minimum point of
the Lasso loss function. The argument \code{tol} must be chosen to transcend the gap
between \eqn{n\lambda} and \eqn{|X^\top_i(y - X\hat\beta)|}. If the tolerance \code{tol}
is chosen too large, some irrelevant features may be included and thus the subsequent facial
reduction procedure would be decelerated; on the other hand, if the tolerance is chosen too 
small, some relevant features may be ignored which leads to a wrong strictly feasible representation. 
In conclusion, a too tight tolerance is more detrimental than a too loose one. An accurate
\eqn{\hat\beta} (\code{beta}) and a slightly loose \code{tol} are always helpful.
}
\examples{
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
y = X \%*\% beta + epsilon
# find one particular solution 
lambda = 0.01
model = glmnet::glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, 
                       intercept = FALSE, standardize = FALSE)
beta = as.numeric( model$beta )

# find the equi-correlation set
equi_index_lasso(X, y, lambda, beta, tol = 1e-5)

}
