\name{netgsaexDAG}
\alias{netgsaexDAG}
\docType{data}
\title{Toy example for using \code{NetGSA} on a two-class test with directed acyclic graphs}
\description{
This dataset contains an example for using Network-based Gene Set Analysis on a two-sample test with directed acyclic graphs.
}
\usage{data("netgsaexDAG")}
\format{
  A list with components
  \describe{
    \item{\code{A}}{A list ot two weighted adjacency matrices.}
    \item{\code{x}}{The \eqn{p \times n}{p x n} data matrix.}
    \item{\code{y}}{The vector of class indicators of length \eqn{n}.}
    \item{\code{B}}{The npath by \eqn{p} indicator matrix for pathways.}
  }
}
\references{
Shojaie, A., & Michailidis, G. (2010). Network enrichment analysis in complex experiments. Statistical applications in genetics and molecular biology, 9(1), Article 22. \url{http://www.ncbi.nlm.nih.gov/pubmed/20597848}.

Shojaie, A., & Michailidis, G. (2009). Analysis of gene sets based on the underlying regulatory network. Journal of Computational Biology, 16(3), 407-426. \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131840/}
}
\examples{
data(netgsaexDAG)
A = netgsaexDAG$A
B = netgsaexDAG$B
x = netgsaexDAG$x
y = netgsaexDAG$y
fitDAG = NetGSA(A, x, y, B, lklMethod="REML", directed=TRUE)
}
\keyword{datasets}
