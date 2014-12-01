\name{covsel}
\alias{covsel}
\title{Covariance selection}
\description{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty.
}
\usage{
covsel(X, zero = NULL, one = NULL, lambda, rho=0.01, verbose = FALSE, eps = 1e-08)
}
\arguments{
  \item{X}{The \eqn{n \times p}{n x p} data matrix.}
  \item{zero}{(Optional) indices of entries of the matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere. The matrix must be symmetric.}
  \item{one}{(Optional) indices of entries of the matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero} and needs to be symmetric.}
  \item{lambda}{(Non-negative) numeric scalar representing the regularization parameter for lasso. }
  \item{rho}{(Non-negative) numeric scalar representing the regularization parameter for estimating the weights in the inverse covariance matrix.}
  \item{verbose}{Whether to print out information as estimation proceeds. Default = \code{FALSE}.}
  \item{eps}{Numeric scalar \eqn{>=0}, indicating the tolerance level for differentiating zero and non-zero edges: entries \eqn{<} \code{eps} will be set to 0. }
}
\details{
The function \code{covsel} performs constrained estimation 
of sparse inverse covariance (concerntration) matrices using a lasso (L1) penalty, as described in Ma, Shojaie and Michailidis (2014). 
Two sets of constraints determine subsets of entries of the inverse covariance matrix that should be exactly zero (the option \code{zero} argument), or should take non-zero values (option \code{one} argument). The remaining entries will be estimated from data.

The arguments \code{one} and/or \code{zero} can come from external knowledge on the 0-1 structure of underlying concerntration matrix, such as a list of edges and/or non-edges learned frm available databases. Then the function \code{edgelist2adj} can be used to first construct \code{one} and/or \code{zero}.

\code{covsel} estimates both the support (0-1 structure) of the concerntration matrix, or equivalently, the adjacency matrix of the corresponding Gaussian graphical model, for a given tuning parameter, \code{lambda}; and the concerntration matrix with diagonal entries set to 0, or equivalently, the weighted adjacency matrix.
The weighted adjacency matrix is estimated using maximum likelihood based on the estimated support. The parameter \code{rho} controls the amount of regularization used in the maximum likelihood step. A small \code{rho} is recommended, as a large value of \code{rho} may result in too much regularization in the maximum likelihood estimation, thus further penalizing the support of the weighted adjacency matrix.  
Note this function is suitable only for estimating the adjacency matrix of a undirected graph.  

This function is closely related to \code{NetGSA}, which requires the weighted adjacency matrix as input. When the user does not have complete information on the weighted adjacency matrix, but has data (\code{X}, not necessarily the same as the \code{x} in \code{NetGSA}) and external information (\code{one} and/or \code{zero}) on the adjacency matrix, then \code{covsel} can be used to estimate the remaining interactions in the adjacency matrix using the data.
Further, when it is anticipated that the adjacency matrices under conditions 1 and 2 are different, and data from both conditions are available, the user needs to run \code{covsel} twice to obtain estimates of the adjacency matrices under each condition.

The algorithm used in \code{covsel} is based on \code{glmnet} and \code{glasso}. Please refer to \code{glmnet} and \code{glasso} for computational details.
}
\value{
A list with components
\item{Adj}{The estimated adjacency matrix of dimension \eqn{p \times p}{p x p}.}
\item{wAdj}{The estimated weighted adjacency matrix of dimension \eqn{p \times p}{p x p}, with diagonal entries set to 0.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted. \url{http://arxiv.org/abs/1411.7919}.
}
\author{
Jing Ma
}
\seealso{
\code{\link{edgelist2adj}}, \code{\link{cv.covsel}}, \code{\link{glmnet}}, \code{\link{glasso}}
}
\examples{
library(MASS)
library(glmnet)
library(glasso)
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

## Covariance selection without external information
fit <- covsel(x, lambda = 0.2)

## Covariance selection with external information
##-Not run-
#oneMat = edgelist2adj(file="edgelist.txt", vertex.names=1:p, mode="undirected")
#zeroMat = edgelist2adj(file="nonedgelist.txt", vertex.names=1:p, mode="undirected")
#fit2 <- covsel(x, zero=zeroMat, one=oneMat, lambda = 0.2)
}
