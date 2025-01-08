#' @keywords internal
#' @noRd
.link_prepare <- function(object, ...) {
  
  # 0) Define local_env as THIS function's environment
  local_env <- environment()
  
  # 1) Retrieve the igraph network(s)
  g <- netw(object)

  # 2) Build the edge data
  if (inherits(object, "mgnet")) {
    edges <- igraph::as_data_frame(g, "edges")
    tbl_from <- edges %>%
      dplyr::select(from) %>%
      dplyr::left_join(gather_taxa(object), by = dplyr::join_by(from == taxa_id)) %>%
      dplyr::rename(taxa_id = from)
    tbl_to <- edges %>%
      dplyr::select(to) %>%
      dplyr::left_join(gather_taxa(object), by = dplyr::join_by(to == taxa_id)) %>%
      dplyr::rename(taxa_id = to)
    
  } else {
    edges <- g %>%
      purrr::map(\(x) igraph::as_data_frame(x, "edges")) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
    tbl_from <- edges %>%
      dplyr::select(mgnet, from) %>%
      dplyr::left_join(gather_taxa(object), by = dplyr::join_by(mgnet, from == taxa_id)) %>%
      dplyr::rename(taxa_id = from)
    tbl_to <- edges %>%
      dplyr::select(mgnet, to) %>%
      dplyr::left_join(gather_taxa(object), by = dplyr::join_by(mgnet, to == taxa_id)) %>%
      dplyr::rename(taxa_id = to)
  }

  # 3) Create auxiliary functions
  # from(var) => pull var from tbl_from
  from <- function(var) {
    var_sym <- rlang::ensym(var)
    dplyr::pull(tbl_from, !!var_sym)
  }
  # to(var) => pull var from tbl_to
  to <- function(var) {
    var_sym <- rlang::ensym(var)
    dplyr::pull(tbl_to, !!var_sym)
  }
  # one(expr) => OR of expr in tbl_from and tbl_to
  either <- function(expr) {
    expr_quo <- rlang::enquo(expr)
    from_vals <- rlang::eval_tidy(expr_quo, data = tbl_from)
    to_vals   <- rlang::eval_tidy(expr_quo, data = tbl_to)
    from_vals | to_vals
  }
  # both(expr) => AND of expr in tbl_from and tbl_to
  both <- function(expr) {
    expr_quo <- rlang::enquo(expr)
    from_vals <- rlang::eval_tidy(expr_quo, data = tbl_from)
    to_vals   <- rlang::eval_tidy(expr_quo, data = tbl_to)
    from_vals & to_vals
  }
  # neither(expr) => NOR (none side)
  neither <- function(expr) {
    expr_quo <- rlang::enquo(expr)
    from_vals <- rlang::eval_tidy(expr_quo, data = tbl_from)
    to_vals   <- rlang::eval_tidy(expr_quo, data = tbl_to)
    !(from_vals | to_vals)
  }
  # one(expr) => XOR, returns TRUE if ONLY one side is TRUE
  one <- function(expr) {
    expr_quo <- rlang::enquo(expr)
    from_vals <- rlang::eval_tidy(expr_quo, data = tbl_from)
    to_vals   <- rlang::eval_tidy(expr_quo, data = tbl_to)
    from_vals != to_vals  
  }
  
  # 4) Capture and re-root user expressions so they see THIS environment
  quos <- rlang::enquos(...)
  quos <- lapply(quos, function(q) {
    rlang::new_quosure(
      expr = rlang::quo_get_expr(q),
      env  = local_env   
    )
  })
  
  # 5) Return a list of all the pieces you might need
  list(
    graph    = g,
    edges    = edges,
    tbl_from = tbl_from,
    tbl_to   = tbl_to,
    from     = from,
    to       = to,
    either   = either,
    both     = both,
    neither  = neither,
    one      = one,
    quos     = quos
  )
}


# SELECT_LINK
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' @title Get Selected Links from an mgnet Object
#'
#' @description
#' Retrieves the numeric indices of links previously selected by [select_link()].
#'
#' @param object An \code{mgnet} or \code{mgnetList} object.
#' @return A numeric vector (or named list of numeric vectors, for \code{mgnetList}) 
#'   containing the indices of selected links.
#' @aliases get_selected_links,mgnet-method get_selected_links,mgnetList-method
#' @export
setGeneric("get_selected_links", function(object) standardGeneric("get_selected_links"))

setMethod("get_selected_links", "mgnet", function(object) {
  return(attr(object, "selected_links"))
})

setMethod("get_selected_links", "mgnetList", function(object) {
  sapply(object, get_selected_links, USE.NAMES = TRUE, simplify = FALSE)
})

#' @title Deselect Links in an mgnet Object
#'
#' @description
#' Clears the \code{selected_links} attribute from an \code{mgnet} or \code{mgnetList} object,
#' effectively removing any current link selection.
#'
#' @param object An \code{mgnet} or \code{mgnetList} object.
#' @return The same object, with no \code{selected_links} attribute.
#' @aliases deselect_link,mgnet-method deselect_link,mgnetList-method
#' @export
setGeneric("deselect_link", function(object) standardGeneric("deselect_link"))

setMethod("deselect_link", "mgnet", function(object) {
  attr(object, "selected_links") <- NULL
  return(object)
})

setMethod("deselect_link", "mgnetList", function(object) {
  for(i in seq_along(object)) object[[i]] <- deselect_link(object[[i]])
  return(object)
})


#' @title Select Links in an mgnet Object Based on Conditions
#'
#' @description
#' This method \strong{highlights} links in the `mgnet` object according to specified
#' conditions and stores the indices of these selected links in the `selected_links`
#' attribute of the object.
#'
#' In particular, you can use `select_link()` \emph{inside} a pipeline that calls
#' [mutate_netw()] to focus on only the selected edges in subsequent network operations.
#'
#' @param object An `mgnet` or `mgnetList` object containing a network (edges).
#' @param ... Conditions to be applied for selecting links, passed as expressions.
#'
#' @details
#' Inside these expressions, you can call the local helper functions documented in
#' [helper_link()], such as \code{from(var)}, \code{to(var)}, or \code{either(expr)}.
#' These helpers allow you to reference and filter based on node-level metadata from
#' the “from” or “to” side of each link.
#'
#' @return 
#' The same `mgnet` (or `mgnetList`) object, updated with a `selected_links`
#' attribute containing the numeric indices of the filtered links.
#'
#' @seealso 
#' [mutate_link()] for modifying link attributes,  
#' [get_selected_links()] and [deselect_link()] for managing selection,  
#' [helper_link()] for details on the local helper functions.
#'
#' @importFrom rlang enquos
#' @importFrom dplyr mutate filter pull
#' @aliases select_link,mgnet-method select_link,mgnetList-method
#' @export
setGeneric("select_link", function(object, ...) standardGeneric("select_link"))

setMethod("select_link", "mgnet", function(object, ...) {
  
  # Ensure the network is available
  if (miss_slot(object, "netw")) stop("Error: No network available.")

  # 1) Internal .link_prepare(), returning a list of everything necessary for links.
  setup <- .link_prepare(object, ...)
  
  # 2) graph, the edges, the helpers, etc.
  g       <- setup$graph
  edges   <- setup$edges
  quos    <- setup$quos
  
  object <- deselect_link(object)

  selected_links <- edges %>%
    dplyr::mutate("_internal_" = row_number()) %>%
    dplyr::filter(!!!quos) %>%
    dplyr::pull("_internal_")

  # Store the filtered links as an attribute
  attr(object, "selected_links") <- selected_links
  return(object)
})

setMethod("select_link", "mgnetList", function(object, ...) {
  
  # Ensure the network is available
  if (miss_slot(object, "netw", "any")) stop("Error: No network available in at least one mgnet element.")
  
  # 1) Internal .link_prepare(), returning a list of everything necessary for links.
  setup <- .link_prepare(object, ...)
  
  # 2) graph, the edges, the helpers, etc.
  g       <- setup$graph
  edges   <- setup$edges
  quos    <- setup$quos
  
  object <- deselect_link(object)
  
  selected_links <- edges %>%
    dplyr::group_by(mgnet) %>%
    dplyr::mutate("_internal_" = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!!quos) %>%
    dplyr::select(mgnet, "_internal_")
  
  selected_links <- split(selected_links, selected_links$mgnet)
  selected_links <- lapply(selected_links, \(x){
    x$mgnet <- NULL
    return(x[["_internal_"]])
  })
  
  # Store the filtered links as an attribute
  for(i in names(object)){
    attr(object[[i]], "selected_links") <- selected_links[[i]]
  }
  return(object)
})


# GROUP LINK
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' @title Retrieve Link Grouping in an mgnet Object
#'
#' @description
#' Retrieves the integer group IDs (or list of group IDs for each sub-object) 
#' previously assigned by [group_link()].
#'
#' @param object An \code{mgnet} or \code{mgnetList} object.
#' @return For an \code{mgnet}, a numeric vector indicating the group ID of each link.
#'   For an \code{mgnetList}, a list of numeric vectors, one per sub-object.
#'
#' @seealso [group_link()] for assigning these group IDs, and [ungroup_link()] 
#'   for removing them.
#'
#' @aliases get_grouped_link,mgnet-method get_grouped_link,mgnetList-method
#' @export
setGeneric("get_grouped_link", function(object) standardGeneric("get_grouped_link"))

#' @rdname get_grouped_link
setMethod("get_grouped_link", "mgnet", function(object) {
  attr(object, "link_groups") 
})

#' @rdname get_grouped_link
setMethod("get_grouped_link", "mgnetList", function(object) {
  lapply(object, get_grouped_link)
})


#' @title Check if Links Are Grouped in an mgnet Object
#'
#' @description
#' Determines whether an `mgnet` or `mgnetList` object currently has grouped links,
#' which would influence behavior inside [mutate_link()].
#'
#' @param object An \code{mgnet} or \code{mgnetList} object.
#' @return A logical value: \code{TRUE} if links are grouped, \code{FALSE} otherwise.
#'         For an `mgnetList`, returns \code{TRUE} only if all sub-objects are grouped.
#' @aliases is_link_grouped,mgnet-method is_link_grouped,mgnetList-method
#' @export
setGeneric("is_link_grouped", function(object) standardGeneric("is_link_grouped"))

setMethod("is_link_grouped", "mgnet", function(object) {
  return(!is.null(get_grouped_link(object)))
})

setMethod("is_link_grouped", "mgnetList", function(object) {
  return(all(sapply(object, is_link_grouped)))
})


#' @title Ungroup Links in an mgnet Object
#'
#' @description
#' Removes any existing link grouping in an `mgnet` or `mgnetList` object, clearing 
#' the \code{link_groups} attribute so that link operations within [mutate_link()] 
#' are no longer grouped.
#'
#' @param object An \code{mgnet} or \code{mgnetList} object.
#' @return The same object, with no \code{link_groups} attribute.
#' @aliases ungroup_link,mgnet-method ungroup_link,mgnetList-method
#' @export
setGeneric("ungroup_link", function(object) standardGeneric("ungroup_link"))

#' @rdname ungroup_link
setMethod("ungroup_link", "mgnet", function(object) {
  attr(object, "link_groups") <- NULL
  object
})

#' @rdname ungroup_link
setMethod("ungroup_link", "mgnetList", function(object) {
  for(i in seq_along(object)) object[[i]] <- ungroup_link(object[[i]])
  object
})


#' @title Group Links in an mgnet Object
#'
#' @description
#' Groups the links (edges) in an `mgnet` or `mgnetList` object based on user-specified
#' expressions, storing an integer group ID for each link in the `link_groups` attribute.
#' These group IDs affect link-level operations within [mutate_link()], allowing you
#' to perform grouped calculations (e.g., grouped summaries or conditional logic) on
#' edge attributes.
#'
#' In particular, you can use `group_link()` \emph{inside} a pipeline with [mutate_link()]
#' to ensure subsequent operations are evaluated by group. 
#'
#' @param object An `mgnet` or `mgnetList` object containing a network (edges).
#' @param ... dplyr-style grouping expressions. For instance, \code{from(family)} or
#'   \code{to(genus)}. See [helper_link()] for details on these local helper functions.
#'
#' @details
#' The grouping expressions can reference node-level metadata in the “from” or “to” side
#' of each link. Once grouped, any subsequent operations in [mutate_link()] can leverage
#' these group IDs to compute grouped statistics or conditionals at the edge level.
#'
#' @return
#' The same `mgnet` (or `mgnetList`) object, updated with a `link_groups` attribute
#' that indicates each link's group ID (or a list of IDs in the case of `mgnetList`).
#'
#' @seealso
#' [mutate_link()] for grouped edge mutations,
#' [ungroup_link()] for clearing group information,
#' [helper_link()] for the local helper functions (e.g., `from()`, `to()`, etc.).
#'
#' @importFrom dplyr group_by group_indices
#' @aliases group_link,mgnet-method group_link,mgnetList-method
#' @export
setGeneric("group_link", function(object, ...) standardGeneric("group_link"))


setMethod("group_link", "mgnet", function(object, ...) {
  
  # Ensure the network is available
  if (miss_slot(object, "netw")) stop("Error: No network available.")
  
  # 1) Prepare edges & expressions via .link_prepare()
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  # 2) Optionally clear any previous grouping
  object <- ungroup_link(object)
  
  # 3) Group the edges using the user expressions
  groups_idx <- edges %>%
    dplyr::group_by(!!!quos) %>%
    dplyr::group_indices()
  
  # 4) Store the grouping info in an attribute
  attr(object, "link_groups") <- groups_idx
  
  # 6) Return the updated object
  object
})


setMethod("group_link", "mgnetList", function(object, ...) {
  
  # Ensure the network is available
  if (miss_slot(object, "netw", "any")) stop("Error: No network available in at least one mgnet element.")
  
  # 1) Prepare with .link_prepare(), which should return edges containing a "mgnet" column
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  # 2) Optionally clear any previous grouping
  object <- ungroup_link(object)
  
  # 3) Group the edges using the user expressions
  groups_idx <- edges %>%
    dplyr::group_by(!!!quos)
  
  groups_idx <- tibble::tibble(
    mgnet = groups_idx$mgnet,
    '_internal_' = dplyr::group_indices(groups_idx)
  )
  
  groups_idx <- split(groups_idx, groups_idx$mgnet)
  groups_idx <- lapply(groups_idx, \(x){
    x$mgnet <- NULL
    return(x[["_internal_"]])
  })
  
  # Store the filtered links as an attribute
  for(i in names(object)){
    attr(object[[i]], "link_groups") <- groups_idx[[i]]
  }
  
  object
})


#' @title Link Helper Functions
#'
#' @description
#' These **helper functions** are designed to manage node-level metadata for the edges
#' (links) in an \code{mgnet} or \code{mgnetList} object. They are primarily intended 
#' for use \emph{inside} link-based methods such as \code{\link{mutate_link}}, 
#' \code{\link{select_link}}, and \code{\link{group_link}}. In those contexts, they
#' have access to the "from" and "to" node metadata for each link.
#'
#' Below is an overview of each helper:
#'
#' \describe{
#'   \item{\code{from(var)}}{
#'     **Retrieves** the metadata column \code{var} from the **first** (or "from") node 
#'     of the link.
#'   }
#'
#'   \item{\code{to(var)}}{
#'     **Retrieves** the metadata column \code{var} from the **second** (or "to") node 
#'     of the link.
#'   }
#'
#'   \item{\code{either(expr)}}{
#'     **Evaluates** \code{expr} in both the "from" and "to" metadata tables for each link, 
#'     returning \code{TRUE} if \emph{either} side satisfies \code{expr}. Conceptually,
#'     it's a logical OR of the two sides.
#'   }
#'
#'   \item{\code{both(expr)}}{
#'     Similar to \code{either(expr)} but returns \code{TRUE} only if \emph{both} sides 
#'     satisfy the condition (logical AND).
#'   }
#'
#'   \item{\code{neither(expr)}}{
#'     Returns \code{TRUE} if \emph{neither} side satisfies \code{expr} (logical NOR). 
#'     In other words, it is \code{TRUE} only if \emph{both} sides evaluate to \code{FALSE}.
#'   }
#'
#'   \item{\code{one(expr)}}{
#'     Returns \code{TRUE} if \emph{exactly one} side satisfies \code{expr}—a logical XOR 
#'     of the "from" side and "to" side. 
#'   }
#' }
#'
#' They \emph{cannot} be called at the top level by name. Instead, they
#' become available \emph{only} within expressions captured by \code{\link{mutate_link}},
#' \code{\link{select_link}}, or \code{\link{group_link}}.
#'
#'
#' @aliases from to either both neither one
#' @export
helper_link <- function() {
  # This function does nothing, 
  # but its Roxygen doc provides a public help page for the link helper functions.
  invisible(NULL)
}

