# graph-coercion.R
# Coercion of mgnet / mgnets objects to tidygraph

#' Coerce mgnet / mgnets objects to a tidygraph tbl_graph
#'
#' These methods convert \code{mgnet} and \code{mgnets} objects into
#' \code{tidygraph::tbl_graph} objects, ready to be used with
#' \code{ggraph::ggraph()}.
#'
#' For a single \code{mgnet} object, the method wraps the internal
#' network (as returned by \code{netw()}) into a \code{tbl_graph}.
#'
#' For a \code{mgnets} object (a collection of \code{mgnet} objects),
#' all internal graphs are combined via \code{igraph::disjoint_union()},
#' after adding a \code{mgnet} column to vertices and edges containing
#' the name (or index) of the originating element. This makes it easy to
#' facet plots by \code{mgnet} in \code{ggraph}.
#'
#' The \code{tidygraph} package is listed in \code{Suggests}; an error is
#' raised at runtime if it is not installed.
#'
#' @param object An object of class \code{mgnet} or \code{mgnets}.
#' @param selected Logical or other selector passed to \code{netw()} to
#'   determine which network (or subset of the network) should be used.
#' @param ... Additional arguments passed to \code{tidygraph::as_tbl_graph()}.
#'
#' @return A \code{tidygraph::tbl_graph} object.
#'
#' @aliases to_tbl_graph,mgnet-method to_tbl_graph,mgnets-method
#' @export
setGeneric("to_tbl_graph", function(object, selected = TRUE, ...) {
  standardGeneric("to_tbl_graph")
})

# -------------------------------------------------------------------------
# to_tbl_graph method for mgnet
# -------------------------------------------------------------------------

#' @rdname to_tbl_graph
#' @export
setMethod(
  "to_tbl_graph", "mgnet",
  function(object, selected = TRUE, ...) {
    
    if (!requireNamespace("tidygraph", quietly = TRUE))
      cli::cli_abort("{.pkg tidygraph} is required for `as_tbl_graph()`.")
    
    if (!requireNamespace("ggraph", quietly = TRUE))
      cli::cli_abort("{.pkg ggraph} is required to visualize the network using {.fn ggraph::ggraph}.")
    
    # Delegate the extraction of the igraph object to netw()
    g <- netw(object, selected = selected)
    
    tidygraph::as_tbl_graph(g, ...)
  }
)

# -------------------------------------------------------------------------
# to_tbl_graph method for mgnets
# -------------------------------------------------------------------------

#' @rdname to_tbl_graph
#' @export
setMethod(
  "to_tbl_graph", "mgnets",
  function(object, selected = TRUE, ...) {
    
    if (!requireNamespace("tidygraph", quietly = TRUE))
      cli::cli_abort("{.pkg tidygraph} is required for `as_tbl_graph()`.")
    
    if (!requireNamespace("ggraph", quietly = TRUE))
      cli::cli_abort("{.pkg ggraph} is required to visualize the network using {.fn ggraph::ggraph}.")
    
    # Handle empty mgnets: return an empty tbl_graph
    if (length(object) == 0L) {
      nodes <- tibble::tibble(.name = character())
      edges <- tibble::tibble()
      return(
        tidygraph::tbl_graph(
          nodes    = nodes,
          edges    = edges,
          directed = FALSE
        )
      )
    }
    
    # Use list names 
    mg_names <- names(object)
    
    # Build one labelled igraph per element
    graph_list <- lapply(seq_along(object), function(i) {
      # Extract the i-th network via netw()
      g <- netw(object[[i]], selected = selected)
      mg_lbl <- mg_names[i]
      
      # Tag vertices with mgnet label
      if (igraph::vcount(g) > 0L) {
        g <- igraph::set_vertex_attr(
          g,
          name  = "mgnet",
          value = rep(mg_lbl, igraph::vcount(g))
        )
      }
      
      # Tag edges with mgnet label
      if (igraph::ecount(g) > 0L) {
        g <- igraph::set_edge_attr(
          g,
          name  = "mgnet",
          value = rep(mg_lbl, igraph::ecount(g))
        )
      }
      
      g
    })
    
    # If there is only one graph, avoid the disjoint_union overhead
    if (length(graph_list) == 1L) {
      g_union <- graph_list[[1L]]
    } else {
      # Duplicate vertex names across graphs are expected here:
      # each component is distinguished by the `mgnet` attribute.
      g_union <- suppressWarnings(
        do.call(igraph::disjoint_union, graph_list)
      )
    }
    
    tidygraph::as_tbl_graph(g_union, ...)
  }
)

