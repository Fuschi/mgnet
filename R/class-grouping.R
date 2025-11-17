#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#==============================================================================#
# Grouping for `mgnet` / `mgnets`
#==============================================================================#

#' Grouping helpers for `mgnet` / `mgnets`
#'
#' Set or get grouping variables on `mgnet` / `mgnets` objects, with semantics
#' similar to `dplyr::group_by()` / `dplyr::ungroup()`.
#'
#' @details
#' The grouping state is stored as an attribute named "mgnet_groups" on the
#' object (for both `mgnet` and `mgnets`).
#'
#' Semantics:
#' * `group_mgnet(x, a, b)` replaces the current groups with `c("a","b")`.
#' * `group_mgnet(x, a, .add = TRUE)` adds `"a"` to the existing groups.
#' * `group_mgnet(x)` with no grouping vars returns `x` unchanged (keep current).
#' * `group_mgnet(x, NULL)` clears the groups (like `group_by(NULL)`).
#' * `ungroup_mgnet(x)` clears all groups.
#' * `ungroup_mgnet(x, a)` removes `"a"` from the current groups (if present).
#' * `is_mgnet_grouped(x)` returns TRUE if grouping vars are stored.
#'
#' Validation:
#' For `mgnet`, valid fields come from `meta_taxa_vars(object)`.
#' For `mgnets`, valid fields come from `meta_taxa_vars(object, .fmt = "unique")`.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param ... One or more variable names (bare or strings). May also include a
#' single `NULL` to drop all groups. In `ungroup_mgnet()`, names in `...` are
#' removed from the current grouping vector.
#' @param .add Logical; if `TRUE`, add to existing groups; if `FALSE`, replace.
#'
#' @return
#' - `group_mgnet()` returns the modified object.
#' - `ungroup_mgnet()` returns the modified object.
#' - `get_group_mgnet()` returns a character vector (possibly length-0).
#' - `is_mgnet_grouped()` returns logical.
#'
#' @name group-mgnet
NULL

# -----------------------------------------------------------------------------#
# Getters
# -----------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("get_group_mgnet", function(object) standardGeneric("get_group_mgnet"))

#' @rdname group-mgnet
#' @export
setMethod("get_group_mgnet", "mgnet", function(object) {return(attr(object, "mgnet_groups"))})

#' @rdname group-mgnet
#' @export
setMethod("get_group_mgnet", "mgnets", function(object) {return(attr(object, "mgnet_groups"))})

#' @importFrom rlang enquos as_name quo_is_null
NULL

# -----------------------------------------------------------------------------#
# Setters (group_mgnet)
# -----------------------------------------------------------------------------#

#' @noRd
.group_mgnet <- function(object, ..., .add = FALSE){
  qs <- rlang::enquos(..., .ignore_empty = "all")
  
  # No arguments or NULL remove all groups
  if (length(qs) == 0L || (length(qs) == 1L && rlang::quo_is_null(qs[[1L]]))) {
    attr(object, "mgnet_groups") <- NULL
    return(object)
  }
  
  # Convert quosures in strings
  vars <- unique(vapply(qs, rlang::as_name, character(1)))
  
  # Validate the arguments
  available <- meta_taxa_vars(object, .fmt = "unique")
  missing_vars <- setdiff(vars, available)
  if (length(missing_vars)) {
    cli::cli_abort(c(
      "x" = "Unknown grouping variable{?s}: {.val {missing_vars}}.",
      "v" = "Available variables are: {.val {available}}."
    ))
  }
  
  if (isTRUE(.add)) {
    current <- get_group_mgnet(object)
    attr(object, "mgnet_groups") <- unique(c(current, vars))
  } else {
    attr(object, "mgnet_groups") <- vars
  }
  object
}

#' @rdname group-mgnet
#' @export
setGeneric("group_mgnet", function(object, ..., .add = FALSE) standardGeneric("group_mgnet"))

#' @rdname group-mgnet
#' @export
setMethod("group_mgnet", signature(object = "mgnet"),
          function(object, ..., .add = FALSE) .group_mgnet(object, ..., .add = .add))

#' @rdname group-mgnet
#' @export
setMethod("group_mgnet", signature(object = "mgnets"),
          function(object, ..., .add = FALSE) .group_mgnet(object, ..., .add = .add))

# -----------------------------------------------------------------------------#
# Ungroup
# -----------------------------------------------------------------------------#

.ungroup_mgnet <- function(object, ...){
  qs <- rlang::enquos(..., .ignore_empty = "all")
  if (length(qs) == 0L) {
    attr(object, "mgnet_groups") <- NULL
    return(object)
  }
  drop <- unique(vapply(qs, rlang::as_name, character(1)))
  keep <- setdiff(get_group_mgnet(object), drop)
  attr(object, "mgnet_groups") <- keep
  object
}

#' @rdname group-mgnet
#' @export
setGeneric("ungroup_mgnet", function(object, ...) standardGeneric("ungroup_mgnet"))

#' @rdname group-mgnet
#' @export
setMethod("ungroup_mgnet", signature(object = "mgnet"), function(object, ...) {
  .ungroup_mgnet(object, ...)
})

#' @rdname group-mgnet
#' @export
setMethod("ungroup_mgnet", signature(object = "mgnets"), function(object, ...) {
  .ungroup_mgnet(object, ...)
})

# -----------------------------------------------------------------------------#
# is_grouped
# -----------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("is_mgnet_grouped", function(object) standardGeneric("is_mgnet_grouped"))

#' @rdname group-mgnet
#' @aliases is_mgnet_grouped,mgnet-method
#' @export
setMethod("is_mgnet_grouped", "mgnet", function(object) {
  !is.null(get_group_mgnet(object))
})

#' @rdname group-mgnet
#' @aliases is_mgnet_grouped,mgnets-method
#' @export
setMethod("is_mgnet_grouped", "mgnets", function(object) {
  all(sapply(object, is_mgnet_grouped))
})
