\name{cv.covsel}
\alias{cv.covsel}
\title{Cross-validation for \code{covsel}}
\description{Performs k-fold cross validation for \code{covsel}.}
\usage{
cv.covsel(X, zero = NULL, one = NULL, lambda, nfolds = NULL, verbose = FALSE)
}
\arguments{
  \item{X}{The \eqn{n \times p}{n x p} data matrix as in \code{covsel}.}
  \item{zero}{(Optional) indices of entries of the matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere. The matrix must be symmetric.}
  \item{one}{(Optional) indices of entries of the matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero} and needs to be symmetric.
  	}
  \item{lambda}{(Non-negative) user-supplied lambda sequence.}
  \item{nfolds}{Number of folds. Default is 5-fold cross validation. }
  \item{verbose}{Whether to print out information as estimation proceeds. Default=\code{FALSE}.}
}
\details{
	The function splits the data \code{X} into \code{nfolds} and runs \code{covsel} \code{nfolds} times to compute the fit with each of the folds omitted. The error is accumulated and the average error over the folds is computed. 
}
\value{
	A list with components
	\item{lambda}{The values of \code{lambda} used.}
	\item{cve}{The mean cross-validated error(s), a vector of \code{length(lambda)}.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted. \url{http://arxiv.org/abs/1411.7919}.
}
\author{
Jing Ma
}
\seealso{
	\code{\link{covsel}}
}
\examples{
library(MASS)
library(glmnet)
library(cvTools)
set.seed(1)

## Generate the covariance matrix for the AR(1) process 
phi <- 0.5
p <- 50
n <- 50
Sigma <- diag(rep(1,p))
Sigma <- phi^(abs(row(Sigma)-col(Sigma)))/(1-phi^2)

## The inverse covariance matrix is sparse
Omega <- solve(Sigma)

## Generate multivariate normal data
x <- mvrnorm(n, mu=rep(0, p), Sigma=Omega)

## Covariance selection with external information
cv.fit <- cv.covsel(x, lambda=seq(0.1,0.3,0.1), nfolds=5)
}
