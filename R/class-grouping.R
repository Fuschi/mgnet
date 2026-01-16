#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#==============================================================================#
# Grouping for `mgnet` / `mgnets`
#==============================================================================#

#' Grouping helpers for `mgnet` / `mgnets`
#'
#' Set, drop, or inspect grouping variables stored on `mgnet` / `mgnets`,
#' with semantics similar to `dplyr::group_by()` / `dplyr::ungroup()`.
#'
#' @details
#' Grouping variables are stored as an attribute named \code{"mgnet_groups"}.
#'
#' Semantics:
#' * `group_mgnet(x, a, b)` replaces groups with `c("a","b")`.
#' * `group_mgnet(x, a, .add = TRUE)` adds `"a"` to existing groups.
#' * `group_mgnet(x)` returns `x` unchanged.
#' * `group_mgnet(x, NULL)` clears all groups.
#' * `ungroup_mgnet(x)` clears all groups.
#' * `ungroup_mgnet(x, a)` removes `"a"` from the current groups.
#' * `get_group_mgnet(x)` returns the grouping vector (or `NULL`).
#' * `is_mgnet_grouped(x)` returns `TRUE` if grouping is active.
#'
#' Validation:
#' Grouping variables must be present in `meta_taxa_vars(object, .fmt = "unique")`.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param ... Grouping variable names (bare or strings). In `group_mgnet()`,
#'   a single `NULL` clears all groups.
#' @param .add Logical; if `TRUE`, add to existing groups; if `FALSE`, replace.
#'
#' @return
#' * `group_mgnet()` / `ungroup_mgnet()` return the modified object.
#' * `get_group_mgnet()` returns a character vector or `NULL`.
#' * `is_mgnet_grouped()` returns a logical scalar.
#'
#' @name group-mgnet
NULL

#' @importFrom rlang enquos as_name quo_is_null
NULL

#------------------------------------------------------------------------------#
# get_group_mgnet
#------------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("get_group_mgnet", function(object) standardGeneric("get_group_mgnet"))

#' @rdname group-mgnet
#' @export
setMethod("get_group_mgnet", "mgnet", function(object) {
  attr(object, "mgnet_groups", exact = TRUE)
})

#' @rdname group-mgnet
#' @export
setMethod("get_group_mgnet", "mgnets", function(object) {
  attr(object, "mgnet_groups", exact = TRUE)
})

#------------------------------------------------------------------------------#
# group_mgnet
#------------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("group_mgnet", function(object, ..., .add = FALSE) standardGeneric("group_mgnet"))

#' @noRd
.group_mgnet <- function(object, ..., .add = FALSE) {
  qs <- rlang::enquos(..., .ignore_empty = "all")
  
  # No args -> keep current
  if (length(qs) == 0L) return(object)
  
  # Single NULL -> clear all groups
  if (length(qs) == 1L && rlang::quo_is_null(qs[[1L]])) {
    attr(object, "mgnet_groups") <- NULL
    return(object)
  }
  
  # Convert quosures to names
  vars <- unique(vapply(qs, rlang::as_name, character(1)))
  
  # Validate names
  available <- meta_taxa_vars(object, .fmt = "unique")
  missing   <- setdiff(vars, available)
  
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Unknown grouping variable{?s}: {.val {missing}}.",
        "v" = "Available variables are: {.val {available}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  # Add or replace
  if (isTRUE(.add)) {
    current <- get_group_mgnet(object)
    vars <- unique(c(if (is.null(current)) character(0) else current, vars))
  }
  
  attr(object, "mgnet_groups") <- vars
  object
}

#' @rdname group-mgnet
#' @export
setMethod("group_mgnet", "mgnet",
          function(object, ..., .add = FALSE) .group_mgnet(object, ..., .add = .add))

#' @rdname group-mgnet
#' @export
setMethod("group_mgnet", "mgnets",
          function(object, ..., .add = FALSE) .group_mgnet(object, ..., .add = .add))

#------------------------------------------------------------------------------#
# ungroup_mgnet
#------------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("ungroup_mgnet", function(object, ...) standardGeneric("ungroup_mgnet"))

#' @noRd
.ungroup_mgnet <- function(object, ...) {
  qs <- rlang::enquos(..., .ignore_empty = "all")
  
  if (length(qs) == 0L) {
    attr(object, "mgnet_groups") <- NULL
    return(object)
  }
  
  drop <- unique(vapply(qs, rlang::as_name, character(1)))
  current <- get_group_mgnet(object)
  if (is.null(current)) return(object)
  
  keep <- setdiff(current, drop)
  attr(object, "mgnet_groups") <- if (length(keep) == 0L) NULL else keep
  object
}

#' @rdname group-mgnet
#' @export
setMethod("ungroup_mgnet", "mgnet", function(object, ...) {
  .ungroup_mgnet(object, ...)
})

#' @rdname group-mgnet
#' @export
setMethod("ungroup_mgnet", "mgnets", function(object, ...) {
  .ungroup_mgnet(object, ...)
})

#------------------------------------------------------------------------------#
# is_mgnet_grouped
#------------------------------------------------------------------------------#

#' @rdname group-mgnet
#' @export
setGeneric("is_mgnet_grouped", function(object) standardGeneric("is_mgnet_grouped"))

#' @rdname group-mgnet
#' @export
setMethod("is_mgnet_grouped", "mgnet", function(object) {
  g <- get_group_mgnet(object)
  is.character(g) && length(g) > 0L
})

#' @rdname group-mgnet
#' @export
setMethod("is_mgnet_grouped", "mgnets", function(object) {
  g <- get_group_mgnet(object)
  is.character(g) && length(g) > 0L
})
