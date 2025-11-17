#------------------------------------------------------------------------------#
#                                 CHECK EXPRESSIONS                            #
#------------------------------------------------------------------------------#

#' Validate a single expression for forbidden dplyr "current-group" calls with omic abundance variables
#'
#' This helper inspects one quosure (or bare expression). It throws an error
#' if it finds both:
#'  - a forbidden dplyr "current-group" function call (e.g., `n()`, `dplyr::n()`)
#'  - and an abundance-related variable (`abun`, `rela`, `norm`) in the same expression.
#'
#' The detection of forbidden functions is based on a regex that matches *calls*
#' (i.e., requires an opening parenthesis), which prevents false positives when
#' a symbol is used as a variable (e.g., `n = ...`).
#'
#' @param quo A quosure from `rlang::enquos(...)` or a bare expression.
#' @param forbidden_functions Character vector of forbidden function names.
#' @param disallowed_variables Character vector of abundance-related variables.
#' @param expr_name Optional label for the expression (used in the error message).
#'
#' @return Invisibly returns `NULL`. Throws with `cli::cli_abort()` on violation.
#' @keywords internal
.validate_abun_expr <- function(
    quo,
    forbidden_functions  = c("n", "cur_group", "cur_group_id", "cur_group_rows", "cur_column"),
    disallowed_variables = c("abun", "rela", "norm"),
    expr_name = NULL
) {
  # Extract underlying expression from a quosure if needed
  expr <- if (inherits(quo, "quosure")) rlang::quo_get_expr(quo) else quo
  
  # Turn the expression into a single string for simple regex checks
  txt <- paste(deparse(expr), collapse = " ")
  
  # Build a regex that matches a *call* to any forbidden function, optionally namespaced.
  # Example matches: "n(", "dplyr::n(", "cur_group(", "pkg.sub::cur_group_rows ("
  fn_pat <- paste0(
    "(?:[[:alnum:]_.]+::)?(",
    paste(forbidden_functions, collapse = "|"),
    ")\\s*\\("
  )
  has_forbidden_call <- grepl(fn_pat, txt, perl = TRUE)
  
  # Gather all symbol names mentioned in the expression (variables, bare names)
  all_names_in_expr <- all.names(expr)
  has_disallowed_name <- length(intersect(disallowed_variables, all_names_in_expr)) > 0
  
  # Also detect .data$abun and .data[['rela']] forms via a lightweight regex
  has_disallowed_dotdata <- grepl(
    "\\.data\\s*\\$\\s*(abun|rela|norm)|\\.data\\s*\\[\\[\\s*['\"](abun|rela|norm)['\"]\\s*\\]\\]",
    txt, perl = TRUE
  )
  
  # If both conditions are met, abort immediately with a clear cli message
  if (has_forbidden_call && (has_disallowed_name || has_disallowed_dotdata)) {
    if (is.null(expr_name) || identical(expr_name, "")) {
      expr_name <- paste0(deparse(expr, width.cutoff = 60), collapse = " ")
    }
    # Build short, informative strings
    # (We don't extract the exact matched function name here to keep things simple.)
    cli::cli_abort(c(
      "Forbidden combination detected.",
      "x" = "Expression {.code {expr_name}} uses a dplyr current-group function (e.g. {.code n()}, {.code cur_group()}) together with omic abundance variables ({.code abun}, {.code rela}, {.code norm}).",
      "i" = "These cannot be used together. Compute any group context first, then reference abundance variables outside that context."
    ))
  }
  
  invisible(NULL)
}

#' Check for Forbidden Functions and Disallowed Variables in Expressions
#'
#' This internal function checks each expression to ensure that forbidden
#' dplyr "current-group" function calls and abundance-related variables do not
#' coexist in the same expression. If they do, an error is thrown immediately.
#'
#' @param expressions A list of expressions captured by `rlang::enquos(...)`.
#'
#' @return Invisibly returns `NULL`. Throws with `cli::cli_abort()` on the first violation.
#' @keywords internal
check_forbidden_expressions <- function(expressions) {
  for (i in seq_along(expressions)) {
    .validate_abun_expr(
      quo       = expressions[[i]],
      expr_name = names(expressions)[i]
    )
  }
  invisible(NULL)
}



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