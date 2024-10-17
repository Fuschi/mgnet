#' Compute Signed Layout for Graphs Considering Only Positive Edges
#'
#' This function generates a graph layout that only considers positive-weighted edges from a given `igraph` or `mgnet` object. 
#' By default, the layout uses the Fruchterman-Reingold algorithm (`with_fr()`), but it can accept other layout functions from `igraph`.
#'
#' @param object An `igraph` or `mgnet` or `mgnetList` object. If the input is an `mgnet` object, the layout is applied to the network (`netw`) slot. If the input is an `mgnetList` the method returns an named list with all the layouts.
#' @param layout A layout specification to be used for computing the graph layout. Defaults to the `Fruchterman-Reingold` layout 
#'        (`with_fr()`), but any layout function from `igraph` can be provided.
#'
#' @details
#' This function focuses on positive-weighted edges in the graph, generating a layout based only on these connections. 
#' If the graph is weighted, edges with positive weights are retained, and negative or zero-weight edges are ignored. 
#' If the graph is unweighted, the function treats all edges as positive.
#'
#' When no layout function is provided, the Fruchterman-Reingold algorithm is used. However, the user can specify 
#' any valid layout from `igraph` to control how the layout is computed.
#'
#' This function is useful for visualizing graphs where positive interactions (e.g., cooperation, mutualism) 
#' are of primary interest.
#'
#' @return A matrix with two columns (for 2D layouts) or three columns (for 3D layouts), containing the coordinates of each vertex.
#'
#' @examples
#' # Create an example igraph object
#' g <- igraph::make_ring(10)
#' E(g)$weight <- c(1, -1, 1, 1, -1, 1, 1, -1, 1, 1)
#'
#' # Compute layout with only positive edges
#' coords <- layout_signed(g)
#'
#' # Plot the graph using the computed layout
#' plot(g, layout = coords)
#'
#' @aliases layout_signed,igraph-method layout_signed,mgnet-method layout_signed,mgnetList-method
#' @importFrom igraph subgraph.edges E is_weighted layout_ with_fr
#' @export
setGeneric("layout_signed", function(object, layout = igraph::with_fr()) standardGeneric("layout_signed"))

setMethod("layout_signed","igraph",function(object, layout = igraph::with_fr()){
  
  if(igraph::is_weighted(object)){
    
    # Check for positive edges
    positive_edges <- which(E(object)$weight > 0)
    object <- igraph::subgraph.edges(object, eids = positive_edges, delete.vertices = FALSE)
    
  } 
  
  # Apply the layout function, defaulting to with_fr() if no layout is provided
  layout_matrix <- igraph::layout_(object, layout)
  
  return(layout_matrix)
})

setMethod("layout_signed","mgnet",function(object, layout = igraph::with_fr()){
  
  layout_signed(netw(object), layout = layout)
  
})

setMethod("layout_signed","mgnetList",function(object, layout = igraph::with_fr()){
  
  sapply(object, \(x) layout_signed(netw(x), layout = layout),
         simplify = FALSE, USE.NAMES = TRUE)
  
})
