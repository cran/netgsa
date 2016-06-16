\name{netgsaex3}
\alias{netgsaex3}
\docType{data}
\title{Toy example for using \code{NetGSA} on a three-class test with undirected networks}
\description{
This dataset contains an example for using Network-based Gene Set Analysis on a multi-sample test with undirected networks.
}
\usage{data("netgsaex3")}
\format{
  A list with components
  \describe{
    \item{\code{A}}{A list ot three weighted adjacency matrices.}
    \item{\code{x}}{The \eqn{p \times n}{p x n} data matrix.}
    \item{\code{y}}{The vector of class indicators of length \eqn{n}.}
    \item{\code{B}}{The npath by \eqn{p} indicator matrix for pathways.}
  }
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information. \url{http://arxiv.org/abs/1411.7919}.
}
\examples{
data(netgsaex3)
A = netgsaex3$A
B = netgsaex3$B
x = netgsaex3$x
y = netgsaex3$y

# -Not-run-
# fit = NetGSA(A, x, y, B, lklMethod="REML", directed=FALSE)
}
\keyword{datasets}
