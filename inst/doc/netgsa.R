## ----setup, include=FALSE-----------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)

## ----data, warning = FALSE, message = FALSE-----------------------------------
library(netgsa)
library(graphite)
library(data.table)
data("breastcancer2012")
ls()

## -----------------------------------------------------------------------------
AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)

## -----------------------------------------------------------------------------
head(rownames(x))

## ----ER status----------------------------------------------------------------
table(group)

## ---- echo = FALSE------------------------------------------------------------
sample_edges <- data.table(base_gene_src = c("7534", "8607"), base_id_src = c("ENTREZID", "ENTREZID"), base_gene_dest = c("8607", "7534"), base_id_dest = c("ENTREZID", "ENTREZID"))
sample_edges

## ----pathways-----------------------------------------------------------------
paths <- graphite::pathways('hsapiens','kegg')
paths[[1]]
head(nodes(paths[[1]]))

## -----------------------------------------------------------------------------
pathways_mat[1:5,7, drop = FALSE]

## ---- echo = FALSE------------------------------------------------------------
as.character(graphite::pathwayDatabases()[graphite::pathwayDatabases()$species == "hsapiens","database"])

## -----------------------------------------------------------------------------
database_search <- obtainEdgeList(rownames(x), c("kegg", "reactome"))
network_info <- prepareAdjMat(x, group, database_search,
                                         cluster = TRUE, file_e = "edgelist.txt", 
                                         file_ne = "nonedgelist.txt")

## -----------------------------------------------------------------------------
network_info[["Adj"]][[1]][[1]][7:9,7:9]

## -----------------------------------------------------------------------------
length(network_info[["Adj"]][[1]])

## -----------------------------------------------------------------------------
pathway_tests_rehe <- NetGSA(network_info[["Adj"]], x, group, pathways_mat, 
                             lklMethod = "REHE", sampling = TRUE, 
                             sample_n = 0.25, sample_p = 0.25)

## ---- eval = FALSE------------------------------------------------------------
#  plot.NetGSA(pathway_tests_rehe)

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
    knitr::include_graphics("cyto_pathway_network_default.png")

## ---- eval = FALSE------------------------------------------------------------
#  formatPathways(pathway_tests_rehe, "ErbB signaling pathway")

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("cyto_erbb.png")

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("cyto_legend.png")

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("igraph_cyto_pathway_network_default.png")

## ---- eval = FALSE------------------------------------------------------------
#  # Format the "Neurotrophin signaling pathway" using the "degree-circle" layout
#  formatPathways(pathway_tests_rehe, "Neurotrophin signaling pathway",
#               graph_layout = "degree-circle")

## ---- eval = FALSE------------------------------------------------------------
#  RCy3::setCurrentNetwork("Pathway Network")
#  edge_weights <- RCy3::getTableColumns(table = "edge", columns = "weight")
#  RCy3::setEdgeLineWidthMapping("weight", c(min(edge_weights), max(edge_weights)), c(1,5))

## ---- eval = FALSE------------------------------------------------------------
#  plot(pathway_tests_rehe)

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("igraph_pathway_network_default.png")

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("igraph_legend.png")

## ---- eval = FALSE------------------------------------------------------------
#  layout_fun <- function(ig) igraph::layout_randomly(ig)
#  plot(pathway_tests_rehe, graph_layout = layout_fun)

## ---- eval = FALSE------------------------------------------------------------
#  zoomPathway(pathway_tests_rehe, "ErbB signaling pathway")

## ---- fig.retina=NULL, out.width=600, echo=FALSE------------------------------
knitr::include_graphics("igraph_erbb.png")

