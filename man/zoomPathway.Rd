\name{zoomPathway}
\alias{zoomPathway}
\title{Zoom in on pathway in igraph}

\description{
Plots the gene to gene interactions for a given pathway in igraph.}

\usage{
zoomPathway(x, pway, graph_layout = NULL)
}

\arguments{
  \item{x}{A NetGSA object returned from calling \code{NetGSA()}}
  \item{pway}{Name of pathway to plot}
  \item{graph_layout}{(Optional) Layout function to pass to igraph plots. This function should only take one parameter (an igraph object). For example one might create a custom layout by setting the spring.length and spring.constant with: \code{my_layout <- function(graph) layout_with_graphopt(graph = graph, spring.length = 1000, spring.constant = 0.00004)}}
}
  
\details{
Generates igraph plot for gene to gene interactions for a given pathway
}

\value{
No return value, called for side effects
}


\references{
Ma, J., Shojaie, A. & Michailidis, G. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics 32(20):165--3174. 
}

\author{
Michael Hellstern
}

\seealso{
\code{\link{plot.NetGSA}}
}

\examples{
\dontrun{
## load the data
data("breastcancer2012_subset")

## consider genes from just 2 pathways
genenames    <- unique(c(pathways[["Adipocytokine signaling pathway"]], 
                         pathways[["Adrenergic signaling in cardiomyocytes"]]))
sx           <- x[match(rownames(x), genenames, nomatch = 0L) > 0L,]

db_edges       <- obtainEdgeList(rownames(sx), databases = c("kegg", "reactome"))
adj_cluster    <- prepareAdjMat(sx, group, databases = db_edges, cluster = TRUE)
out_cluster    <- NetGSA(adj_cluster[["Adj"]], sx, group, 
                         pathways_mat[c(1,2), rownames(sx)], lklMethod = "REHE", sampling = FALSE)
plot(out_cluster)
my_layout <- function(graph) layout_with_graphopt(graph = graph, 
                                                  spring.length = 1000, 
                                                  spring.constant = 0.00004)
zoomPathway(out_cluster, "Adipocytokine signaling pathway", my_layout)
}
}
