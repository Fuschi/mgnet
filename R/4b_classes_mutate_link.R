#' Mutate Link Data for an mgnet object
#'
#' @description
#' `mutate_link,mgnet-method` provides a data-frame-like interface for mutating 
#' edge-level attributes in an igraph network associated with an `mgnet` object. 
#'
#' It can utilize columns from the "from" and "to" nodes via local helpers:
#'
#' - `from(var)`: Pull the column `var` from the node on the "from" side.
#' - `to(var)`: Pull the column `var` from the node on the "to" side.
#' - `either(expr)`: For each edge, returns `TRUE` if \emph{either} side meets `expr`.
#' - `both(expr)`: For each edge, returns `TRUE` if \emph{both} sides meet `expr`.
#' - `neither(expr)`: For each edge, returns `TRUE` if \emph{neither} side meets `expr`.
#' - `one(expr)`: Returns \code{TRUE} if \emph{exactly one} side satisfies `expr`.
#'
#' After constructing (or modifying) columns, the new columns are attached as
#' edge attributes in the igraph object. The `mgnet` object is then returned
#' with an updated network.
#'
#' @param object An `mgnet` object containing an igraph network (retrieved by
#'   `netw(object)`) and taxon data (retrieved by `gather_taxa(object)`).
#' @param ... `dplyr`-style mutate expressions. Inside these expressions, you 
#'   can call the local helpers described above.
#'
#' @return The original `mgnet` object, updated with new or modified edge attributes.
#'
#' @rdname mutate_link
#' @aliases mutate_link,mgnet-method
#' @export
setGeneric("mutate_link", function(object, ...) {standardGeneric("mutate_link")})

setMethod("mutate_link", "mgnet", function(object, ...) {
  # 1) Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g     <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if(is_link_grouped(object)){
    edges <-  edges %>%
      dplyr::mutate(`_internal_` = get_grouped_link(object)) %>%
      dplyr::group_by(`_internal_`)
  }
  
  # 3) Perform the mutate on (potentially) grouped edges
  edges <- dplyr::mutate(edges, !!!quos)
  
  # 4) Attach any new columns as edge attributes
  edge_cols <- setdiff(names(edges), c("from", "to", "_internal_"))
  for (col_name in edge_cols) {
    g <- igraph::set_edge_attr(
      g,
      name  = col_name,
      index = igraph::E(g),
      value = edges[[col_name]]
    )
  }
  
  # 5) Update the igraph object in 'object' and return
  netw(object) <- g
  object
})


setMethod("mutate_link", "mgnetList", function(object, ...) {
  # 1) Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if(is_link_grouped(object)){
    edges <-  edges %>%
      dplyr::mutate(`_internal_` = unlist(get_grouped_link(object), use.names = F)) %>%
      dplyr::group_by(`_internal_`)
  }
  
  # 3) Perform the mutate on (potentially) grouped edges
  edges <- dplyr::mutate(edges, !!!quos)
  edges_names <- names(edges)

  edges <- split(edges, edges$mgnet)
  edges <- sapply(edges, \(x){
    x$mgnet <- NULL
    return(x)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # 4) Attach any new columns as edge attributes
  edge_cols <- setdiff(edges_names, c("from", "to", "_internal_"))
  for (nm in names(object)){
    for (col_name in edge_cols) {
      g[[nm]] <- igraph::set_edge_attr(
        g[[nm]],
        name  = col_name,
        index = igraph::E(g[[nm]]),
        value = edges[[nm]][[col_name]]
      )
    }
    # 5) Update the igraph object in 'object' and return
    netw(object[[nm]]) <- g[[nm]]
  }
  
  object
})

