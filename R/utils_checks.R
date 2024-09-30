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
  reserved_keywords <- c("sample_id", "taxa_id", "comm_id", "abun", "rela", "norm", "mgnet")
  
  # Search reserved keywords in the expression names
  found_keywords <- names(expressions)[names(expressions) %in% reserved_keywords]
  found_keywords <- unique(found_keywords)
  
  # If any reserved keywords were found, raise an error
  if (length(found_keywords) > 0) {
    stop(sprintf("Error: The following reserved keywords can't be modified: %s",
                 paste(found_keywords, collapse = ", ")))
  }
}

#------------------------------------------------------------------------------#
#' Validate Required Variables in Expressions
#'
#' This internal function checks that all variables required by the expressions are present
#' in the object (either in abundance data or metadata). It returns the needed abundance
#' and metadata keys.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param keys_required A character vector of variable names required by the expressions.
#' @param sample_or_taxa A character string, either "sample" or "taxa", indicating the context.
#'
#' @return A list containing:
#'   \describe{
#'     \item{needed_abundance_keys}{Character vector of needed abundance keys (e.g., "abun", "rela", "norm").}
#'     \item{needed_noabundance_keys}{Character vector of needed metadata keys.}
#'   }
#' @keywords internal
validate_required_variables <- function(object, keys_required, sample_or_taxa) {
  
  keys_abundance <- c("abun", "rela", "norm")
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Get available metadata variables
  if (sample_or_taxa == "sample") {
    metadata_columns <- unique(unlist(meta_vars(object)))
  } else {
    metadata_columns <- unique(unlist(taxa_vars(object)))
  }
  
  # Check for missing abundance data
  missing_abundances <- c()
  for (key in keys_abundance) {
    if (key %in% keys_required && length(methods::slot(object, key)) == 0) {
      missing_abundances <- c(missing_abundances, key)
    }
  }
  
  # Check for missing metadata variables
  missing_noabundances <- c()
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundances <- setdiff(needed_noabundance_keys, metadata_columns)
  }
  
  # Generate error message if any keys are missing
  if (length(missing_abundances) > 0 || length(missing_noabundances) > 0) {
    error_message <- "Error: The following required variables in the expressions are missing:\n"
    
    # Add missing abundance-related keys
    if (length(missing_abundances) > 0) {
      error_message <- paste0(error_message, "- abundance-related: ", 
                              paste(missing_abundances, collapse = ", "), "\n")
    }
    
    # Add missing metadata keys
    if (length(missing_noabundances) > 0) {
      error_message <- paste0(error_message, "- ", sample_or_taxa, " metadata: ", 
                              paste(missing_noabundances, collapse = ", "), "\n")
    }
    
    stop(error_message)
  }
  
  return(list(
    needed_abundance_keys = needed_abundance_keys,
    needed_noabundance_keys = needed_noabundance_keys
  ))
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
