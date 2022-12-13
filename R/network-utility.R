################################################################################
################################################################################
# SIGNED LAYOUT
################################################################################
################################################################################
#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.
#' 
#' @param graph network belong to igraph class
#' @param seed random seed for reprodubility
#'
#' @importFrom igraph layout.fruchterman.reingold subgraph.edges E is.igraph
#' @export
layout_signed <- function(graph, seed=123){
  
  if(is.igraph(graph)) stop("graph must belong to igraph class")
  if(!is.numeric(seed)) stop("seed must be numeric")
  
  graph.sub <- subgraph.edges(graph=graph,
                              eids=which(E(graph)$weight>0),
                              delete.vertices=FALSE)
  
  set.seed(seed)
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
}
################################################################################
################################################################################
# END SIGNED LAYOUT
################################################################################
################################################################################




################################################################################
################################################################################
# COMMUNITIES COLORMAP
################################################################################
################################################################################
#' Create color palette for communities
#' 
#' @description User wrapper for distinctColorPalette function of randomcoloR
#' package. Create n+1 distinct color to associate with communities ID. The 
#' first color labeled as 0 is always white for isolated nodes.
#' 
#' @param n positive integer indicated the number of distinct color.
#' @param seed random seed for reproducibility.
#' 
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom grDevices rgb
#' @export
colormap_communities <- function(n=20, seed=123){
  
  if(round(n)!=n | !is.numeric(n) | n<=0) stop("n must be an integer positive number")
  
  set.seed(1)
  colormap <- randomcoloR::distinctColorPalette(k=n)
  colormap <- c(grDevices::rgb(1,1,1,.8),colormap)
  names(colormap) <- as.character(0:n)
  
  return(colormap)
}
################################################################################
################################################################################
# END COMMUNITIES COLORMAP
################################################################################
################################################################################




################################################################################
################################################################################
# TAXONOMY COLORMAP
################################################################################
################################################################################
#' Create color palette for taxonomy
#' 
#' @description User wrapper for distinctColorPalette function of randomcoloR
#' package. Create named vector with taxaID as name and colors as values. If are
#' present duplicated in the vector the function choice only unique values.
#' 
#' @param taxaID character vector with all taxonomic classification.
#' @param seed random seed for reproducibility.
#' 
#' @importFrom randomcoloR distinctColorPalette
#' @export
colormap_taxonomy <- function(taxaID, seed=123){
  
  if(!is.character(taxaID)) stop("taxaID must be character")
  if(!is.null(dim(taxaID))) stop("taxaID can't be a matrix")
  if(!is.numeric(seed)) stop("seed must be numeric")

  taxaID <- unique(taxaID)
  n <- length(taxaID)
  
  set.seed(seed)
  colormap <- randomcoloR::distinctColorPalette(k=n)
  names(colormap) <- taxaID
  
  return(colormap)
}