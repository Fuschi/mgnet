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
  
  # Define the reserved keywords
  reserved_keywords <- c("sample_id", "taxa_id", "comm_id", "abun", "rela", "norm", "mgnet", "meta", "taxa", ".")
  
  # Search reserved keywords in the expression names
  found_keywords <- names(expressions)[names(expressions) %in% reserved_keywords]
  found_keywords <- unique(found_keywords)
  
  # If any reserved keywords were found, raise an error
  if (length(found_keywords) > 0) {
    stop(sprintf("Error: The following reserved keywords can't be modified: %s",
                 paste(found_keywords, collapse = ", ")))
  }
}

#' Check for Forbidden Functions and Disallowed Variables in Expressions
#'
#' This internal function checks each expression to ensure that forbidden functions
#' and disallowed variables do not coexist in any expression. If they do, an error is thrown.
#'
#' @param expressions A list of expressions captured by `rlang::enquos(...)`.
#'
#' @return Throws an error if any expression contains both a forbidden function and a disallowed variable; otherwise, nothing is returned.
#' @keywords internal
check_forbidden_expressions <- function(expressions) {
  
  # Forbidden functions and disallowed variables
  forbidden_functions <- c("n", "cur_group", "cur_group_id", "cur_group_rows", "cur_column")
  disallowed_variables <- c("abun", "rela", "norm")
  
  check_expression_detailed <- function(expr) {
    
    # Extract all names from the expression (functions and variables)
    all_names_in_expr <- all.names(expr)
    
    # Check if any forbidden functions are present
    found_forbidden_function <- intersect(forbidden_functions, all_names_in_expr)
    
    # Check if any disallowed variables are present
    found_disallowed_variable <- intersect(disallowed_variables, all_names_in_expr)
    
    # If both a forbidden function and a disallowed variable are present, return details
    if (length(found_forbidden_function) > 0 && length(found_disallowed_variable) > 0) {
      return(list(
        has_violation = TRUE,
        forbidden_function = found_forbidden_function,
        disallowed_variable = found_disallowed_variable
      ))
    }
    
    # No violation
    return(list(has_violation = FALSE))
  }
  
  # Check each expression for violations and return detailed errors
  for (i in seq_along(expressions)) {
    expr_name <- names(expressions)[i]
    expr_content <- rlang::quo_get_expr(expressions[[i]])
    
    # Get the violation details
    check_result <- check_expression_detailed(expr_content)
    
    # If there is a violation, stop and show detailed error message
    if (check_result$has_violation) {
      stop(sprintf(
        "Error: Expression '%s' contains a current group function from `dplyr` ('%s') and an abundance-related variable from `mgnet` ('%s'). These cannot be used together. See the `dplyr` context for more information.",
        expr_name,
        paste(check_result$forbidden_function, collapse = ", "),
        paste(check_result$disallowed_variable, collapse = ", ")
      ))
    }
  }
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
