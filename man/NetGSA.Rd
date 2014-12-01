\name{NetGSA}
\alias{NetGSA}
\title{Network-based Gene Set Analysis}
\description{
Tests the significance of pre-defined sets of genes (pathways) with respect to an outcome variable, such as the condition indicator (e.g. cancer vs. normal, etc.), based on the underlying biological network.
}
\usage{
NetGSA(A1, A2, x, y, B, lklMethod = c("REML", "ML"), directed = FALSE, eta = 0.1,
       lim4kappa = 500)
}
\arguments{
  \item{A1}{The weighted adjacency matrix for condition 1.}
  \item{A2}{The weighted adjacency matrix for condition 2.}
  \item{x}{The \eqn{p \times n}{p x n} data matrix.}
  \item{y}{Vector of class indicators of length \eqn{n}.}
  \item{B}{The npath by \eqn{p} indicator matrix for pathways.}
  \item{lklMethod}{Method used for likelihood calculation: options are \code{ML} (maximum likelihood) or \code{REML} (restricted maximum likelihood).}
  \item{directed}{Whether the networks are directed.}
  \item{eta}{Approximation limit for the Influence matrix. See 'Details'.}
  \item{lim4kappa}{Limit for condition number (used to adjust \code{eta}). See 'Details'.}
}
\details{
The function \code{NetGSA} carries out a Network-based Gene Set Analysis, using the method described in Shojaie and Michailidis (2009) and Shojaie and Michailidis (2010). It differs from Gene Set Analysis (Efron and Tibshirani, 2007) in that it incorporates the underlying biological networks.

The NetGSA method is formulated in terms of a mixed linear model. Let \eqn{X} represent the rearrangement of data \code{x} into an \eqn{np \times 1}{np x 1} column vector.  
\deqn{X=\Psi \beta + \Pi \gamma + \epsilon}
where \eqn{\beta} is the vector of fixed effects, \eqn{\gamma} and \eqn{\epsilon} are random effects and random errors, respectively. The underlying biological networks are encoded in the weighted adjacency matrices \code{A1} and \code{A2}, which determine the
influence matrix under each condition. The influence matrices further determine the design matrices \eqn{\Psi} and \eqn{\Pi} in the mixed linear model. Formally, the influence matrix under each condition represents the effect of each gene on all the other genes in the network and
is calculated from the adjacency matrix (\code{A1} or \code{A2}). A small value of \code{eta} is used to make sure that the influence matrices are well-conditioned (i.e. their condition numbers are bounded by \code{lim4kappa}.)

The problem is then to test the null hypothesis \eqn{\ell\beta = 0}{l\beta = 0} vs. the alterernative \eqn{\ell\beta \neq 0}{l\beta != 0}, where \eqn{\ell}{l} is a contrast vector, optimally defined through the underlying networks. 
The test statistic \eqn{T} for each gene set is then a function of \eqn{\beta}, variances of \eqn{\gamma} and \eqn{\epsilon}, the constrast vector \eqn{\ell}{l} and the underlying biological network(s) in both conditions. Under the null hypothesis, \eqn{T} has approximately a \eqn{t}-distribution, whose degrees of freedom are estimated using the Satterthwaite approximation method. 
The fixed effects \eqn{\beta} are estimated by generalized least squares, and the estimate depends on estimates of the variance components of \eqn{\gamma} and \eqn{\epsilon}. The variance components (\eqn{\sigma^2_{\epsilon}} and \eqn{\sigma^2_{\gamma}}) are estimated using Newton's method based on the profiling out \eqn{\sigma_{\epsilon}}.

This function can deal with both directed and undirected networks, which are specified via the option \code{directed}. Note \code{NetGSA} uses slightly different procedures to calculate the influence matrices for directed and undirected networks. 
In the case of undirected networks, the user can still apply \code{NetGSA} if only partial information on the adjacency matrices is available. The function \code{covsel} provides one way to estimate the weighted adjacency matrices from data based on available network information. 
}
\value{
A list with components
\item{beta}{Vector of fixed effects of length \eqn{2p}, of which the first half is for condition 1 and the second half for conditin 2.}
\item{df}{Degrees of freedom for the test statistics.}
\item{p.value}{P-values for gene sets (pathways).}
\item{s2.epsilon}{Variance of the random errors \eqn{\epsilon}.}
\item{s2.gamma}{Variance of the random effects \eqn{\gamma}.}
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted. \url{http://arxiv.org/abs/1411.7919}

Shojaie, A., & Michailidis, G. (2010). Network enrichment analysis in complex experiments. Statistical applications in genetics and molecular biology, 9(1), Article 22.

Shojaie, A., & Michailidis, G. (2009). Analysis of gene sets based on the underlying regulatory network. Journal of Computational Biology, 16(3), 407-426.
}
\author{
Ali Shojaie and Jing Ma
}
\seealso{
\code{\link{edgelist2adj}}, \code{\link{covsel}}
}
\examples{
set.seed(1)
library(igraph)
data(netgsaex)

A1 = netgsaex$A1
A2 = netgsaex$A2
B = netgsaex$B
x = netgsaex$x
y = netgsaex$y

##Visualize the networks
par(mar = c(0.5, 0.5, 3, 0.5))
plot(netgsaex$g.alt, vertex.size = 5, vertex.label = NA, main="Network - alt")

par(mar = c(0.5, 0.5, 3, 0.5))
plot(netgsaex$g.null, vertex.size = 5, vertex.label = NA, main="Network - null")

out = NetGSA(A1, A2, x, y, B, lklMethod = "REML")
out2 = NetGSA(A1, A2, x, y, B, lklMethod = "ML")
}
