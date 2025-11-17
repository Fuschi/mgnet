#'@include class-mgnet.R class-mgnets.R class-links.R
NULL

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
#'   `netw(object)`) and taxon data.
#' @param ... `dplyr`-style mutate expressions. Inside these expressions, you 
#'   can call the local helpers described above.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the link grouping
#'        from the object at the end.
#' @param .deselect Logical, default `FALSE`. If `TRUE`, remove the link selection
#'        from the object at the end.
#'
#' @return The original `mgnet` object, updated with new or modified edge attributes.
#'
#' @rdname mutate_link
#' @aliases mutate_link,mgnet-method mutate_link,mgnets-method
#' @importFrom dplyr mutate group_by
#' @importFrom igraph set_edge_attr E
#' @importFrom rlang .data
#' @export
setGeneric("mutate_link", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {standardGeneric("mutate_link")})

setMethod("mutate_link", "mgnet", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  # Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g     <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if (is_link_grouped(object)) {
    edges <- edges %>%
      dplyr::mutate(`_internal_` = get_grouped_link(object)) %>%
      dplyr::group_by(.data[["_internal_"]])
  }
  
  link(object) <- dplyr::mutate(edges, !!!quos) %>% dplyr::ungroup()
  
  if(isTRUE(.ungroup)) object <- ungroup_link(object)
  if(isTRUE(.deselect)) object <- deselect_link(object)
  object
})


setMethod("mutate_link", "mgnets", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  # Prepare the edges, expressions, and local helpers via .link_prepare()
  setup <- .link_prepare(object, ...)
  g <- setup$graph
  edges <- setup$edges
  quos  <- setup$quos
  
  if(is_link_grouped(object)) {
    link_tbl_groups <- object %>%
      purrr::imap(\(x, nm){
      tibble::tibble(mgnet = nm, link_id = link_id(x),
                     `_internal_` = get_grouped_link(x))}) %>%
      purrr::list_rbind()
    
    edges <- edges %>%
      dplyr::left_join(link_tbl_groups, by = c("mgnet", "link_id")) %>%
      dplyr::group_by(.data[["_internal_"]])
  }

  # Perform the mutate on (potentially) grouped edges
  link(object) <- dplyr::mutate(edges, !!!quos) %>% dplyr::ungroup() %>%
    dplyr::select(-tidyselect::any_of("_internal_"))
  
  if(isTRUE(.ungroup)) object <- ungroup_link(object)
  if(isTRUE(.deselect)) object <- deselect_link(object)
  object
})

  