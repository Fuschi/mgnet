#------------------------------------------------------------------------------#
#' Create `link_id` values for a graph
#'
#' @description
#' Generate a character vector of edge identifiers in the form
#' `"<from>--<to>--<index>"`.
#'
#' These identifiers are intended to be unique within the current graph and are
#' used internally to support persistent edge selection.
#'
#' @param g An `igraph` object.
#'
#' @return
#' A character vector of length equal to the number of edges in `g`.
#'
#' @details
#' The generated IDs depend on the current edge order returned by
#' `igraph::as_data_frame(g, what = "edges")`. They are therefore stable within
#' the current graph instance, but may change if the graph is rebuilt or edges
#' are reordered.
#'
#' @keywords internal
#------------------------------------------------------------------------------#
.make_link_id <- function(g) {
  
  # No edges -> return empty character vector
  m <- igraph::ecount(g)
  if (m == 0L) return(character(0L))
  
  # Extract edges as a data frame
  df <- igraph::as_data_frame(g, what = "edges")
  
  # Generate IDs
  ids <- paste(df$from, df$to, seq_len(nrow(df)), sep = "--")
  
  ids
}


#------------------------------------------------------------------------------#
#' Ensure that a graph has a `link_id` edge attribute
#'
#' @description
#' Check whether an `igraph` object already has an edge attribute named
#' `link_id`. If not, generate one with [.make_link_id()] and attach it to the
#' graph.
#'
#' @param g An `igraph` object.
#'
#' @return
#' The input graph, with a guaranteed `link_id` edge attribute whenever edges
#' are present.
#'
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
  
  g
}