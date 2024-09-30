# SIGNED LAYOUT
#------------------------------------------------------------------------------#
#' Signed-Weighted Graph Layout
#'
#' Generates layout coordinates for a graph by considering only positive-weighted edges.
#' This function is applicable to graphs of the `igraph` and `mgnet` classes. It emphasizes 
#' the structure formed by positive interactions by ignoring edges with non-positive weights.
#'
#' @param object A network object belonging to either the `igraph` or `mgnet` class.
#'            The graph must be weighted and contain at least one positive-weighted edge.
#'
#' @return A matrix of vertex coordinates determined by the Fruchterman-Reingold layout algorithm, 
#'         considering only positive-weighted edges.
#'
#' @export
#' @importFrom igraph layout.fruchterman.reingold subgraph.edges E is_weighted
#' @aliases layout_signed,igraph-method layout_signed,mgnet-method
setGeneric("layout_signed", function(object) standardGeneric("layout_signed"))

setMethod("layout_signed","igraph",function(object){
  
  if(!is_weighted(object)) stop("object must be a weighted graph")
  graph.sub <- igraph::subgraph.edges(graph=object,
                                      eids=which(E(object)$weight>0),
                                      delete.vertices=FALSE)
  
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
})

setMethod("layout_signed","mgnet",function(object){
  
  graph.sub <- subgraph.edges(graph=netw(object),
                              eids=which(E(netw(object))$weight>0),
                              delete.vertices=FALSE)
  
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
})

