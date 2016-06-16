\name{netEst.undir}
\alias{netEst.undir}
\title{
Constrained estimation of undirected networks
}
\description{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty.
}
\usage{
netEst.undir(X, zero = NULL, one = NULL, lambda, rho = NULL, weight = NULL, 
             eta = 0, verbose = FALSE, eps = 1e-08)
}
\arguments{
  \item{X}{The \eqn{n \times p}{n x p} data matrix.}
  \item{zero}{(Optional) indices of entries of the matrix to be constrained to be zero. The input should be a matrix of \eqn{p \times p}{p x p}, with 1 at entries to be constrained to be zero and 0 elsewhere. The matrix must be symmetric.}
  \item{one}{(Optional) indices of entries of the matrix to be kept regardless of the regularization parameter for lasso. The input is similar to that of \code{zero} and needs to be symmetric.}
  \item{lambda}{(Non-negative) numeric scalar representing the regularization parameter for lasso. This algorithm only accepts one \code{lambda} at a time.}
  \item{rho}{(Non-negative) numeric scalar representing the regularization parameter for estimating the weights in the inverse covariance matrix.}
  \item{weight}{(Optional) whether to add penalty to known edges. If NULL (default), then the known edges are assumed to be true. If nonzero, then a penalty equal to lambda * weight is added to penalize the known edges to account for possible uncertainty. Only non-negative values are accepted for the weight parameter.  }
  \item{eta}{(Non-negative) a small constant added to the diagonal of the empirical covariance matrix of \code{X} to ensure it is well conditioned. By default, \code{eta} is set to 0.}
  \item{verbose}{Whether to print out information as estimation proceeds. Default = \code{FALSE}.}
  \item{eps}{(Non-negative) numeric scalar indicating the tolerance level for differentiating zero and non-zero edges: entries with magnitude \eqn{<} \code{eps} will be set to 0. }
}
\details{
The function \code{netEst.undir} performs constrained estimation of sparse inverse covariance (concerntration) matrices using a lasso (L1) penalty, as described in Ma, Shojaie and Michailidis (2014). Two sets of constraints determine subsets of entries of the inverse covariance matrix that should be exactly zero (the option \code{zero} argument), or should take non-zero values (option \code{one} argument). The remaining entries will be estimated from data.

The arguments \code{one} and/or \code{zero} can come from external knowledge on the 0-1 structure of underlying concerntration matrix, such as a list of edges and/or non-edges learned frm available databases. Then the function \code{edgelist2adj} can be used to first construct \code{one} and/or \code{zero}.

\code{netEst.undir} estimates both the support (0-1 structure) of the concerntration matrix, or equivalently, the adjacency matrix of the corresponding Gaussian graphical model, for a given tuning parameter, \code{lambda}; and the concerntration matrix with diagonal entries set to 0, or equivalently, the weighted adjacency matrix.
The weighted adjacency matrix is estimated using maximum likelihood based on the estimated support. The parameter \code{rho} controls the amount of regularization used in the maximum likelihood step. A small \code{rho} is recommended, as a large value of \code{rho} may result in too much regularization in the maximum likelihood estimation, thus further penalizing the support of the weighted adjacency matrix.  
Note this function is suitable only for estimating the adjacency matrix of a undirected graph. The \code{weight} parameter allows one to specify whether to penalize the known edges. If known edges obtained from external information contain uncertainty such that some of them are spurious, then it is recommended to use a small positive \code{weight} parameter to select the most probable edges from the collection of known ones. 

This function is closely related to \code{NetGSA}, which requires the weighted adjacency matrix as input. When the user does not have complete information on the weighted adjacency matrix, but has data (\code{X}, not necessarily the same as the \code{x} in \code{NetGSA}) and external information (\code{one} and/or \code{zero}) on the adjacency matrix, then \code{netEst.undir} can be used to estimate the remaining interactions in the adjacency matrix using the data.
Further, when it is anticipated that the adjacency matrices under different conditions are different, and data from different conditions are available, the user needs to run \code{netEst.undir} separately to obtain estimates of the adjacency matrices under each condition.

The algorithm used in \code{netEst.undir} is based on \code{glmnet} and \code{glasso}. Please refer to \code{glmnet} and \code{glasso} for computational details.
}
\value{
A list with components
\item{Adj}{The weighted adjacency matrix (partial correlations) of dimension \eqn{p \times p}{p x p}, with diagonal entries set to 0. This is the matrix that will be used in \code{NetGSA}.}
\item{invcov}{The estimated inverse covariance matrix of dimension \eqn{p \times p}{p x p}.}
\item{lambda}{The values of tuning parameters used.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted. \url{http://arxiv.org/abs/1411.7919}.
}
\author{
Jing Ma
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{edgelist2adj}},\code{\link{bic.netEst.undir}}, \code{\link{glmnet}}, \code{\link{glasso}}
}
\examples{
library(MASS)
library(glmnet)
library(glasso)
set.seed(1)

## Generate the covariance matrix for the AR(1) process 
rho <- 0.5
p <- 100
n <- 100
Sigma <- diag(rep(1,p))
Sigma <- rho^(abs(row(Sigma)-col(Sigma)))/(1-rho^2)

## The inverse covariance matrix is sparse
Omega <- solve(Sigma)

## Generate multivariate normal data
X <- mvrnorm(n, mu=rep(0, p), Sigma=Omega)

## Estimate the network without external information
fit <- netEst.undir(X, lambda = 0.2)

## Estimate the network with external information
##-Not run-
ones = edgelist2adj(file="edgelist.txt", vertex.names=paste0("gene", 1:p), 
mode="undirected")
zeros = edgelist2adj(file="nonedgelist.txt", vertex.names=paste0("gene", 1:p), 
mode="undirected")

fit.undir <- netEst.undir(X, zero=zeros, one=ones, lambda = 0.5)

## Estimate the network when the known edges are not entirely reliable. 
ones[1,10:11] = 1
ones[10:11,1] = 1
fit.undir <- netEst.undir(X, zero=zeros, one=ones, weight = 0.1, lambda = 0.2)
}
