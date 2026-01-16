#------------------------------------------------------------------------------#
#' Check for Reserved Keywords in Expression Names
#'
#' This internal function checks if any reserved keywords are used in the names of expressions.
#' Reserved keywords cannot be modified and should not be used as names in the expressions.
#'
#' @param expressions A list of expressions captured by `rlang::enquos(...)`.
#'
#' @return Throws an error if any reserved keywords are found in the expression names; otherwise, nothing is returned.
#' @keywords internal
check_reserved_keywords <- function(expressions) {
  
  # Normalize names (handles NULL) and drop empty names
  expr_names <- rlang::names2(expressions)
  expr_names <- expr_names[expr_names != ""]
  
  # Intersect with reserved keywords
  found <- unique(expr_names[expr_names %in% .MGNET_RESERVED_COLNAMES])
  
  if (length(found)) {
    cli::cli_abort(c(
      "Reserved keywords cannot be modified.",
      "x" = "Found: {.val {found}}"
    ))
  }
  
  invisible(NULL)
}


#------------------------------------------------------------------------------#
#' Extract Abundance-Related Keys from Expressions
#'
#' This function analyzes expressions to identify any references to specific abundance-related keys within an `mgnet` object. 
#' It checks for mentions of "abun", "rela", or "norm", which represent different types of abundance data.
#'
#' @param expressions A list of expressions, typically created using `rlang::enquos()` in a user-facing function.
#'        These expressions are meant to manipulate or assess the abundance data stored in an `mgnet` object.
#'
#' @return A character vector containing unique abundance-related keys found within the provided expressions.
#'        These keys are essential for further processing or filtering operations in metagenomic data analysis.
#'
#' @details The function parses each expression to extract variables, which are then checked against known abundance-related
#'          keys. It's primarily used internally within functions that need to process or manipulate abundance data based on
#'          user-defined criteria.
#'
#' @keywords internal
#'
#' @importFrom purrr map
#' @importFrom rlang quo_get_expr
get_abundance_keys <- function(expressions){
  
  keys_abuns <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique() %>%
    intersect(c("abun", "rela", "norm"))
  
  return(keys_abuns)
  
}