\name{edgelist2adj}
\alias{edgelist2adj}
\title{Construct an adjacency matrix from an edge list}
\description{
Read the edge list of a graph from a file and construct the adjacency matrix.
}
\usage{
edgelist2adj(file, vertex.names, mode = c("directed", "undirected"))
}
\arguments{
  \item{file}{The connection to read from. This should be a .txt file, with one edge in a line, the two vertices separated by a delimiter.}
  \item{vertex.names}{The names of all vertices in the graph.}
  \item{mode}{Whether the graph to read is directed.}
}
\details{
The function \code{edgelist2adj} accepts a list of edges and constructs the 0-1 adjacency matrix corresponding to the graph. If \code{file} contains an incomplete list of edges, the function determines the actual size of the graph through the vector \code{vertex.names}. 

Although named as \code{edgelist2adj}, the user can also construct a 0-1 matrix corresponding to non-edges, i.e. connections that do not exist between any two vertices.

When \code{file} contains incomplete information, the returned 0-1 adjacency matrix (for edges or non-edges) can be used as input in \code{covsel} to estimate the complete (and weighted) adjacency matrix. 
}
\value{
A 0-1 adjacency matrix of dimension \eqn{p \times p} corresponding to the edges defined in \code{file}, where \eqn{p} is the length of the vector \code{vertex.names}.
}
\references{
Ma, J., Shojaie, A. & Michailidis, G. (2014). Network-based pathway enrichment analysis with incomplete network information, submitted. \url{http://arxiv.org/abs/1411.7919}.
}
\author{
Ali Shojaie and Jing Ma
}
\seealso{
\code{\link{netEst.dir}},\code{\link{netEst.undir}}
}
\examples{
#Read the data
data(edgelist)
data(nonedgelist)

#Generate the .txt files
write.table(edgelist, file="edgelist.txt", row.names=FALSE)
write.table(nonedgelist, file="nonedgelist.txt", row.names=FALSE)

#Read the edge/nonedge list from files
oneMat = edgelist2adj(file="edgelist.txt", vertex.names=paste0("gene", 1:100), 
         mode="undirected")
zeroMat = edgelist2adj(file="nonedgelist.txt", vertex.names=paste0("gene", 1:100), 
          mode="undirected")
}
