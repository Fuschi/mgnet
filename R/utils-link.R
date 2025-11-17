#------------------------------------------------------------------------------#
#' Create link_id values for a graph
#'
#' Generates a character vector of link IDs in the form:
#' "<from>--<to>--<index>".
#'
#' @param g An igraph object.
#'
#' @return A character vector of length equal to the number of edges.
#' @keywords internal
#------------------------------------------------------------------------------#
.make_link_id <- function(g) {
  
  # No edges → return empty character vector
  m <- igraph::ecount(g)
  if (m == 0L) return(character(0L))
  
  # Extract edges as a data.frame
  df <- igraph::as_data_frame(g, what = "edges")
  
  # Generate the IDs
  ids <- paste(df$from, df$to, seq_len(nrow(df)), sep = "--")
  
  return(ids)
}


#------------------------------------------------------------------------------#
#' Ensure the graph has a link_id edge attribute
#'
#' This function checks whether the "link_id" attribute is defined. If missing, 
#' it generates link IDs using `.make_link_id()` and assigns them.
#'
#' @param g An `igraph` object.
#'
#' @return The `igraph` object, with guaranteed link_id in the network if edges exist.
#' @keywords internal
#------------------------------------------------------------------------------#
.ensure_link_id <- function(g) {
  
  # No edges: nothing to do
  if (igraph::ecount(g) == 0L) return(g)
  
  # If link_id already exists, nothing to do
  if (!is.null(igraph::edge_attr(g, "link_id"))) return(g)
  
  # Otherwise generate link_id
  ids <- .make_link_id(g)
  
  # Assign link_id
  g <- igraph::set_edge_attr(
    graph = g,
    name  = "link_id",
    value = ids
  )
  
  return(g)
}


