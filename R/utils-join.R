#' Left-join after pruning overlapping columns in `x` or `y`
#'
#' Removes overlapping non-key columns from either `x` or `y` (controlled by
#' `prefer`) and then performs `dplyr::left_join()`. This avoids `.x`/`.y`
#' suffixes and lets you choose which table's overlapping columns "win".
#'
#' If `prefer = "x"`, columns overlapping between `x` and `y` (excluding the
#' join keys in `x`) are removed from `x` so the result takes those columns
#' from `y`. If `prefer = "y"`, the overlapping non-key columns are removed
#' from `y`, keeping `x`'s versions in the result.
#'
#' @param x A data frame (left table).
#' @param y A data frame (right table).
#' @param by A character vector or named character vector of join keys
#'   (as in `dplyr::left_join()`), e.g. `by = "id"` or `by = c(id_x = "id_y")`.
#' @param prefer One of `"x"` or `"y"`. Choose which table to prune before joining.
#' @param ... Passed to `dplyr::left_join()`.
#' @param verbose Logical; if `TRUE`, prints which columns are pruned.
#'
#' @return A data frame resulting from the left join.
#' @keywords internal
#' @noRd
left_join_prune <- function(x, y, by, prefer = c("x", "y"), ..., verbose = TRUE) {
  prefer <- match.arg(prefer)
  
  # Support by = c("id") or by = c(id_in_x = "id_in_y")
  if (is.null(names(by))) {
    key_x <- by
    key_y <- by
  } else {
    key_x <- names(by)
    key_y <- unname(by)
  }
  
  if (prefer == "x") {
    # Remove from x any non-key columns that also exist in y
    overlap <- intersect(setdiff(names(x), key_x), names(y))
    if (length(overlap)) {
      if (isTRUE(verbose)) {
        message("Pruning from `x`: ", paste(overlap, collapse = ", "))
      }
      x <- dplyr::select(x, -dplyr::all_of(overlap))
    }
  } else { # prefer == "y"
    # Remove from y any non-key columns that also exist in x
    overlap <- intersect(setdiff(names(y), key_y), names(x))
    if (length(overlap)) {
      if (isTRUE(verbose)) {
        message("Pruning from `y`: ", paste(overlap, collapse = ", "))
      }
      y <- dplyr::select(y, -dplyr::all_of(overlap))
    }
  }
  
  dplyr::left_join(x, y, by = by, ...)
}

