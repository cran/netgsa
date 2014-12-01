\name{netgsaex}
\alias{netgsaex}
\docType{data}
\title{Toy example for using \code{NetGSA} on a two-class test}
\description{
This dataset contains an example for using Network-based Gene Set Analysis on a two-sample test.
}
\usage{data(netgsaex)}
\format{
  A list with components
  \describe{
    \item{\code{A1}}{The weighted adjacency matrix for condition 1.}
    \item{\code{A2}}{The weighted adjacency matrix for condition 2.}
    \item{\code{x}}{The \eqn{p \times n}{p x n} data matrix.}
    \item{\code{y}}{The vector of class indicators of length \eqn{n}.}
    \item{\code{B}}{The npath by \eqn{p} indicator matrix for pathways.}
    \item{\code{g.null}, \code{g.alt}}{The networks under the null and alternative hypothesis, respectively. Node colors in the two graphs represent different mean values. Note the network structures between \code{g.null} and \code{g.alt} differ slightly.}
  }
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted.
}
\examples{
data(netgsaex)
}
\keyword{datasets}
