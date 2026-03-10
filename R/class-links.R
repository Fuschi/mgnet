#'@include class-mgnet.R class-mgnets.R

#------------------------------------------------------------------------------#
#' Prepare link-level data and helper functions
#'
#' @description
#' Internal helper used by link-oriented methods such as [select_link()],
#' [group_link()], and [mutate_link()]. It builds an edge-level table and the
#' corresponding `"from"` / `"to"` node-level metadata tables, then exposes a
#' set of local helper functions (`from()`, `to()`, `either()`, `both()`,
#' `neither()`, `one()`) inside captured expressions.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param . Captured expressions passed from a caller such as [select_link()] or
#'   [group_link()].
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `graph`: the underlying graph object, or list of graphs for `mgnets`;
#'   \item `edges`: the edge-level tibble used for downstream filtering/grouping;
#'   \item `tbl_from`: metadata table aligned to the `from` endpoint of each edge;
#'   \item `tbl_to`: metadata table aligned to the `to` endpoint of each edge;
#'   \item `from`, `to`, `either`, `both`, `neither`, `one`: local helper functions;
#'   \item `quos`: captured expressions re-rooted to the local environment.
#' }
#'
#' @keywords internal
#' @noRd
.link_prepare <- function(object, ...) {
  
  # 0) Define local_env as THIS function's environment
  local_env <- environment()
  
  # 1) Retrieve the igraph network(s)
  g <- netw(object, selected = FALSE)
  
  # 2) Build the edge data
  if (inherits(object, "mgnet")) {
    
    edges <- igraph::as_data_frame(g, "edges")
    tbl_from <- edges %>%
      dplyr::select(from) %>%
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = dplyr::join_by(from == taxa_id)) %>%
      dplyr::rename(taxa_id = from)
    tbl_to <- edges %>%
      dplyr::select(to) %>%
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = dplyr::join_by(to == taxa_id)) %>%
      dplyr::rename(taxa_id = to)
    
  } else {
    
    edges <- g %>%
      purrr::map(\(x) igraph::as_data_frame(x, "edges")) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
    tbl_from <- edges %>%
      dplyr::select(mgnet, from) %>%
      dplyr::left_join(taxa(object, .collapse = TRUE), by = dplyr::join_by(mgnet, from == taxa_id)) %>%
      dplyr::rename(taxa_id = from)
    tbl_to <- edges %>%
      dplyr::select(mgnet, to) %>%
      dplyr::left_join(taxa(object, .collapse = TRUE), by = dplyr::join_by(mgnet, to == taxa_id)) %>%
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


#------------------------------------------------------------------------------#
# LINK HELPER DOCUMENTATION
#------------------------------------------------------------------------------#

#' @title Link helper functions
#'
#' @description
#' These helper functions are available inside link-oriented methods such as
#' [mutate_link()], [select_link()], and [group_link()]. They expose node-level
#' metadata aligned to the two endpoints of each edge.
#'
#' They are intended for use with both `mgnet` and `mgnets` objects.
#'
#' @details
#' The following helpers are available inside captured expressions:
#'
#' \describe{
#'   \item{\code{from(var)}}{
#'     Return the values of metadata column `var` for the `"from"` endpoint of
#'     each edge.
#'   }
#'
#'   \item{\code{to(var)}}{
#'     Return the values of metadata column `var` for the `"to"` endpoint of
#'     each edge.
#'   }
#'
#'   \item{\code{either(expr)}}{
#'     Evaluate `expr` on both endpoints and return `TRUE` whenever at least one
#'     side is `TRUE`.
#'   }
#'
#'   \item{\code{both(expr)}}{
#'     Evaluate `expr` on both endpoints and return `TRUE` only when both sides
#'     are `TRUE`.
#'   }
#'
#'   \item{\code{neither(expr)}}{
#'     Evaluate `expr` on both endpoints and return `TRUE` only when neither
#'     side is `TRUE`.
#'   }
#'
#'   \item{\code{one(expr)}}{
#'     Evaluate `expr` on both endpoints and return `TRUE` only when exactly one
#'     side is `TRUE`.
#'   }
#' }
#'
#' These helpers are only available inside expressions captured by
#' [mutate_link()], [select_link()], and [group_link()]. They are not intended
#' to be called directly at top level.
#'
#' Missing values are not coerced to `FALSE`: helper expressions follow standard
#' R logical semantics. As a consequence, if an expression returns `NA`, it will
#' usually be dropped by downstream `dplyr::filter()` calls.
#'
#' @aliases from to either both neither one
#' @export
helper_link <- function() {
  invisible(NULL)
}


#------------------------------------------------------------------------------#
# SELECTED LINKS
#------------------------------------------------------------------------------#

#' @title Retrieve selected link IDs
#'
#' @description
#' Return the values currently stored in the `selected_links` attribute.
#'
#' Selected links are stored as edge `link_id` values, not as positional edge
#' indices.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' For an `mgnet`, an atomic vector of selected `link_id` values, or `NULL` if
#' no selection is currently set.
#'
#' For an `mgnets`, a named list with one element per contained `mgnet`, each
#' element being either a vector of selected `link_id` values, `character(0)`,
#' or `NULL`.
#'
#' @aliases get_selected_links,mgnet-method get_selected_links,mgnets-method
#' @export
setGeneric("get_selected_links", function(object) standardGeneric("get_selected_links"))

#' @rdname get_selected_links
setMethod("get_selected_links", "mgnet", function(object) {
  attr(object, "selected_links")
})

#' @rdname get_selected_links
setMethod("get_selected_links", "mgnets", function(object) {
  sapply(object, get_selected_links, USE.NAMES = TRUE, simplify = FALSE)
})


#' @title Check whether links are selected
#'
#' @description
#' Determine whether an `mgnet` or `mgnets` object currently carries an active
#' link selection.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' For an `mgnet`, `TRUE` if the `selected_links` attribute is set, `FALSE`
#' otherwise.
#'
#' For an `mgnets`, `TRUE` only if **all** contained `mgnet` objects currently
#' carry a link selection. This implements an all-or-none collection-level
#' semantics.
#'
#' @aliases are_selected_links,mgnet-method are_selected_links,mgnets-method
#' @export
setGeneric("are_selected_links", function(object) standardGeneric("are_selected_links"))

#' @rdname are_selected_links
setMethod("are_selected_links", "mgnet", function(object) {
  !is.null(attr(object, "selected_links"))
})

#' @rdname are_selected_links
setMethod("are_selected_links", "mgnets", function(object) {
  all(sapply(object, are_selected_links, USE.NAMES = TRUE, simplify = TRUE))
})


#' @title Clear link selection
#'
#' @description
#' Remove the `selected_links` attribute from an `mgnet` or `mgnets` object.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' The same object, with no active link selection.
#'
#' @aliases deselect_link,mgnet-method deselect_link,mgnets-method
#' @export
setGeneric("deselect_link", function(object) standardGeneric("deselect_link"))

#' @rdname deselect_link
setMethod("deselect_link", "mgnet", function(object) {
  attr(object, "selected_links") <- NULL
  object
})

#' @rdname deselect_link
setMethod("deselect_link", "mgnets", function(object) {
  for (i in seq_along(object)) object[[i]] <- deselect_link(object[[i]])
  object
})


#' @title Select links by condition
#'
#' @description
#' Select edges in an `mgnet` or `mgnets` object according to one or more
#' filtering expressions and store the resulting edge `link_id` values in the
#' `selected_links` attribute.
#'
#' This selection can then be used by downstream methods such as [mutate_link()]
#' or by helpers that explicitly query selected edges.
#'
#' @param object An `mgnet` or `mgnets` object containing a network.
#' @param ... Filtering expressions evaluated on edge-level data. These
#'   expressions may use the local helpers documented in [helper_link()], such
#'   as `from()`, `to()`, `either()`, `both()`, `neither()`, and `one()`.
#'
#' @details
#' For `mgnets`, selection is evaluated on the combined edge table across the
#' collection, but the resulting `link_id` vectors are stored back into each
#' contained `mgnet`.
#'
#' If an individual `mgnet` in the collection has a network but no edges match
#' the selection criteria, its `selected_links` attribute is set to
#' `character(0)`. By contrast, `NULL` means that no selection has been set at
#' all.
#'
#' @return
#' The same object, updated with a `selected_links` attribute containing selected
#' edge `link_id` values.
#'
#' @seealso
#' [mutate_link()] for link-level transformations,
#' [get_selected_links()] and [deselect_link()] for selection management,
#' [helper_link()] for the local helper functions available inside expressions.
#'
#' @aliases select_link,mgnet-method select_link,mgnets-method
#' @importFrom rlang enquos
#' @importFrom dplyr filter pull group_by summarise select all_of
#' @export
setGeneric("select_link", function(object, ...) standardGeneric("select_link"))

#' @rdname select_link
setMethod("select_link", "mgnet", function(object, ...) {
  
  if (miss_netw(object)) {
    cli::cli_abort("No network available for {.cls mgnet} object.")
  }
  
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  object <- deselect_link(object)
  
  selected_links <- edges %>%
    dplyr::filter(!!!quos) %>%
    dplyr::pull("link_id")
  
  attr(object, "selected_links") <- selected_links
  object
})

#' @rdname select_link
setMethod("select_link", "mgnets", function(object, ...) {
  
  if (miss_netw(object, "any")) {
    cli::cli_abort(
      "No network available in at least one element of the {.cls mgnets} object."
    )
  }
  
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  object <- deselect_link(object)
  
  selected_links <- edges %>%
    dplyr::group_by(.data$mgnet) %>%
    dplyr::filter(!!!quos) %>%
    dplyr::select(dplyr::all_of(c("mgnet", "link_id"))) %>%
    dplyr::summarise(link_id = list(link_id), .groups = "drop") %>%
    tibble::deframe()
  
  for (nm in names(object)) {
    g_i <- object[[nm]]@netw
    
    attr(object[[nm]], "selected_links") <-
      if (igraph::ecount(g_i) == 0L) {
        # No edges in the network: store an empty selection (character(0)), not NULL
        character(0)
      } else if (nm %in% names(selected_links)) {
        # Matching edges found: store their link_id values
        selected_links[[nm]]
      } else {
        # Edges exist but none satisfy the filter: empty selection (character(0))
        character(0)
      }
  }
  
  object
})


#------------------------------------------------------------------------------#
# GROUPED LINKS
#------------------------------------------------------------------------------#

#' @title Retrieve link grouping
#'
#' @description
#' Return the values currently stored in the `link_groups` attribute.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' For an `mgnet`, an integer vector giving the group ID of each edge, or `NULL`
#' if no grouping is currently set.
#'
#' For an `mgnets`, a named list of integer vectors, one per contained `mgnet`.
#'
#' @seealso
#' [group_link()] for assigning link groups and [ungroup_link()] for removing
#' them.
#'
#' @aliases get_grouped_link,mgnet-method get_grouped_link,mgnets-method
#' @export
setGeneric("get_grouped_link", function(object) standardGeneric("get_grouped_link"))

#' @rdname get_grouped_link
setMethod("get_grouped_link", "mgnet", function(object) {
  attr(object, "link_groups")
})

#' @rdname get_grouped_link
setMethod("get_grouped_link", "mgnets", function(object) {
  lapply(object, get_grouped_link)
})


#' @title Check whether links are grouped
#'
#' @description
#' Determine whether an `mgnet` or `mgnets` object currently carries link
#' grouping information.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' For an `mgnet`, `TRUE` if the `link_groups` attribute is set, `FALSE`
#' otherwise.
#'
#' For an `mgnets`, `TRUE` only if **all** contained `mgnet` objects currently
#' carry link grouping information. This implements an all-or-none
#' collection-level semantics.
#'
#' @aliases is_link_grouped,mgnet-method is_link_grouped,mgnets-method
#' @export
setGeneric("is_link_grouped", function(object) standardGeneric("is_link_grouped"))

#' @rdname is_link_grouped
setMethod("is_link_grouped", "mgnet", function(object) {
  !is.null(get_grouped_link(object))
})

#' @rdname is_link_grouped
setMethod("is_link_grouped", "mgnets", function(object) {
  all(sapply(object, is_link_grouped, USE.NAMES = TRUE, simplify = TRUE))
})


#' @title Clear link grouping
#'
#' @description
#' Remove the `link_groups` attribute from an `mgnet` or `mgnets` object.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' The same object, with no active link grouping.
#'
#' @aliases ungroup_link,mgnet-method ungroup_link,mgnets-method
#' @export
setGeneric("ungroup_link", function(object) standardGeneric("ungroup_link"))

#' @rdname ungroup_link
setMethod("ungroup_link", "mgnet", function(object) {
  attr(object, "link_groups") <- NULL
  object
})

#' @rdname ungroup_link
setMethod("ungroup_link", "mgnets", function(object) {
  for (i in seq_along(object)) object[[i]] <- ungroup_link(object[[i]])
  object
})


#' @title Group links by condition
#'
#' @description
#' Group edges in an `mgnet` or `mgnets` object according to one or more
#' grouping expressions and store the resulting integer group IDs in the
#' `link_groups` attribute.
#'
#' These group IDs can then be used by downstream methods such as [mutate_link()]
#' to perform grouped edge-level calculations.
#'
#' @param object An `mgnet` or `mgnets` object containing a network.
#' @param ... Grouping expressions evaluated on edge-level data. These
#'   expressions may use the local helpers documented in [helper_link()], such
#'   as `from()`, `to()`, `either()`, `both()`, `neither()`, and `one()`.
#'
#' @details
#' For a single `mgnet`, group IDs are computed on that object's edge table.
#'
#' For `mgnets`, grouping is computed at the **collection level** on the
#' combined edge table. This means that edges from different contained `mgnet`
#' objects can receive the same group ID if they share the same grouping key.
#'
#' After grouping is computed on the combined table, each vector of group IDs is
#' stored back into the corresponding contained `mgnet`.
#'
#' @return
#' The same object, updated with a `link_groups` attribute.
#'
#' @seealso
#' [mutate_link()] for grouped link-level transformations,
#' [ungroup_link()] for clearing group information,
#' [helper_link()] for the local helper functions available inside expressions.
#'
#' @aliases group_link,mgnet-method group_link,mgnets-method
#' @importFrom dplyr group_by group_indices
#' @export
setGeneric("group_link", function(object, ...) standardGeneric("group_link"))

#' @rdname group_link
setMethod("group_link", "mgnet", function(object, ...) {
  
  if (miss_netw(object)) {
    cli::cli_abort("No network available for {.cls mgnet} object.")
  }
  
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  object <- ungroup_link(object)
  
  groups_idx <- edges %>%
    dplyr::group_by(!!!quos) %>%
    dplyr::group_indices()
  
  attr(object, "link_groups") <- groups_idx
  object
})

#' @rdname group_link
setMethod("group_link", "mgnets", function(object, ...) {
  
  if (miss_netw(object, "any")) {
    cli::cli_abort(
      "No network available in at least one element of the {.cls mgnets} object."
    )
  }
  
  setup <- .link_prepare(object, ...)
  edges <- setup$edges
  quos  <- setup$quos
  
  object <- ungroup_link(object)
  
  groups_tbl <- edges %>%
    dplyr::group_by(!!!quos)
  
  groups_tbl <- tibble::tibble(
    mgnet      = groups_tbl$mgnet,
    `_group_id` = dplyr::group_indices(groups_tbl)
  )
  
  groups_idx <- split(groups_tbl, groups_tbl$mgnet)
  groups_idx <- lapply(groups_idx, function(x) {
    x$mgnet <- NULL
    as.integer(x[["_group_id"]])
  })
  
  for (nm in names(object)) {
    g_i <- object[[nm]]@netw
    
    attr(object[[nm]], "link_groups") <-
      if (igraph::ecount(g_i) == 0L) {
        # No edges in the network: store an empty grouping (integer(0))
        integer(0)
      } else if (nm %in% names(groups_idx)) {
        # Group IDs computed for this mgnet
        groups_idx[[nm]]
      } else {
        # Edges exist but no group indices were generated (safety fallback)
        integer(0)
      }
  }
  
  object
})