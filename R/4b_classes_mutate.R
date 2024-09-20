#' Modify and Augment `mgnet` Objects by Transforming `sample` Slot
#'
#' This function dynamically manipulates the `sample` data slot within `mgnet` objects,
#' applying user-defined transformations. It harnesses the comprehensive functionality of the
#' `tidyverse`, particularly `dplyr`, to facilitate these data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object from which data will be modified.
#' This function specifically targets the `sample` slot, which stores metadata for each sample.
#' @param ... Dynamic expressions or functions to be applied to the data, utilizing
#' \code{\link[dplyr]{mutate}} directly. This integration supports all standard and custom
#' `tidyverse` manipulation techniques, ensuring versatile and complex data transformations.
#' @param .by Optional; a character vector indicating the columns used for grouping data
#' before applying transformations. For `mgnet` objects, defaults to "sample_id".
#' For `mgnetList` objects, defaults to \code{c("mgnet", "sample_id")}, where "mgnet"
#' refers to individual `mgnet` objects within the list, ensuring transformations are
#' applied contextually within each sub-object. The use of "taxa_id" as a grouping variable
#' is explicitly prohibited to maintain consistency and prevent data misinterpretation.
#'
#' @details This function is designed to enhance and modify the metadata of samples stored within
#' the `sample` slot using conditional logic and statistical summaries tailored to the needs
#' of the analysis. It is built to integrate seamlessly with the `tidyverse` methods, promoting
#' familiar and powerful data manipulation capabilities directly within this package. Users are
#' encouraged to apply any `tidyverse`-compatible expressions or `dplyr` verbs to achieve the desired
#' data transformation, including but not limited to conditional operations, summarizations, or the
#' use of `dplyr::across`.
#'
#' @return Returns the `mgnet` or `mgnetList` object with updated `sample` slots reflecting
#' the applied transformations. All other structures within the object remain unchanged, ensuring
#' that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_sample,mgnet-method mutate_sample,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate
#' @importFrom tidyr expand_grid any_of
#' @importFrom rlang enquos quo_get_expr syms exprs
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_sample", function(object, ..., .by) {standardGeneric("mutate_sample")})

setMethod("mutate_sample", "mgnet", function(object, ..., .by = "sample_id") {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Ensure there are samples to process
  if (nsample(object) == 0) {
    stop("Error: No samples available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::quos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "sample_id"
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Retrieve sample columns
  sample_columns <- sample_vars(object)
  
  # Check additional elements in .by against sample columns
  additional_by_elements <- setdiff(.by, "sample_id")
  if (length(additional_by_elements) > 0 && any(!additional_by_elements %in% sample_columns)) {
    stop(sprintf("Error: All elements in '.by' except 'sample_id' must be columns in 'sample'. Missing: `%s`",
                 paste(setdiff(additional_by_elements, sample_columns), collapse = "`, `")))
  }
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()

  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")

  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)

  # Prepare to check for missing keys
  missing_keys <- list()

  # Check each abundance key for presence
  for (key in keys_abundance) {
    if (key %in% keys_required && length(methods::slot(object, key)) == 0) {
      missing_keys <- c(missing_keys, key)
    }
  }
  
  # Check for missing non-abundance keys in sample
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundance_keys <- setdiff(needed_noabundance_keys, sample_columns)
    if (length(missing_noabundance_keys) > 0) {
      missing_keys <- c(missing_keys, missing_noabundance_keys)
    }
  }
  
  # VALIDATION CHECK: Ensure no forbidden functions and disallowed variables coexist
  #----------------------------------------------------------------------------#
  
  # Forbidden functions and disallowed variables
  forbidden_functions <- c("n", "cur_group", "cur_group_id", "cur_group_rows", "cur_column")
  disallowed_variables <- c("abun", "rela", "norm")
  
  # Function to check for forbidden functions and disallowed variables in an expression
  check_expression_detailed <- function(expr, forbidden_functions, disallowed_variables) {
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
    check_result <- check_expression_detailed(expr_content, forbidden_functions, disallowed_variables)
    
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
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  if(length(needed_abundance_keys) != 0){
    
    sample_abun_data <- tidyr::expand_grid(sample_id = sample_id(object),
                                           taxa_id = taxa_id(object))
    
    for(abundance_key in needed_abundance_keys){
      sample_abun_data <- sample_abun_data %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id, 
                                names_to = "taxa_id", 
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
    
    if(length(sample(object)) != 0){
      sample_abun_data <- sample_abun_data %>%
        dplyr::left_join(sample(object, .fmt = "tbl"),
                         by = "sample_id")
    }
    
  }
  # END LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  # APPLY MUTATE
  #----------------------------------------------------------------------------#
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(rlang::quo_get_expr(expr))
    
    #expr_name <- names(expressions)[i]
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      result <- sample_abun_data %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        #dplyr::mutate(!!rlang::sym(expr_name) := !!expr) %>%
        dplyr::mutate(!!expr) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("taxa_id", keys_abundance))) %>%
        dplyr::distinct() %>%
        dplyr::arrange(match(sample_id, sample_id(object)))
      
    } else {
      
      # Process expressions not involving abundance-related variables
      result <- sample(object, "tbl") %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        #dplyr::mutate(!!rlang::sym(expr_name) := !!expr) %>%
        dplyr::mutate(!!expr) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(match(sample_id, sample_id(object)))
    }
    
    expr_name <- tail(names(result), 1)
    # Update sample in the object
    if (!expr_name %in% names(sample(object))) {
      # Add new column if it doesn't exist
      sample(object) <- sample(object) %>%
        dplyr::mutate(!!rlang::sym(expr_name) := result[[expr_name]])
    } else {
      # Replace existing column
      sample(object)[[expr_name]] <- result[[expr_name]]
    }
    
  }
  
  validObject(object)
  return(object)
  
})


setMethod("mutate_sample", "mgnetList", function(object, ..., .by = c("mgnet", "sample_id")) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Ensure there are samples to process
  if (any(nsample(object) == 0)) {
    stop("Error: No samples available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::quos(...)
  
  # # Validate that all expressions are named
  # expression_names <- names(expressions)
  # 
  # if (is.null(expression_names) || any(nzchar(expression_names) == FALSE)) {
  #   stop("All expressions must be named. Please provide names for all transformations.")
  # }
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","sample_id")
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Retrieve sample columns
  sample_columns <- unique(unname(unlist(sample_vars(object))))
  
  # Check additional elements in .by against sample columns
  additional_by_elements <- setdiff(.by, c("mgnet","sample_id"))
  if (length(additional_by_elements) > 0 && any(!additional_by_elements %in% sample_columns)) {
    stop(sprintf("Error: All elements in '.by' except 'sample_id' must be columns in 'sample'. Missing: `%s`",
                 paste(setdiff(additional_by_elements, sample_columns), collapse = "`, `")))
  }
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")
  
  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Prepare to check for missing keys
  missing_keys <- list()
  
  # Check each abundance key for presence
  for (i in 1:length(object)){
    mgnet_obj <- object[[i]]
    for (key in keys_abundance) {
      if (key %in% keys_required && length(methods::slot(mgnet_obj, key)) == 0) {
        missing_keys <- c(missing_keys, key)
      }
    }
  }
  
  # Check for missing non-abundance keys in sample
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundance_keys <- setdiff(needed_noabundance_keys, sample_columns)
    if (length(missing_noabundance_keys) > 0) {
      missing_keys <- c(missing_keys, missing_noabundance_keys)
    }
  }
  
  # Report all missing keys, if any
  if (length(missing_keys) > 0) {
    stop(sprintf("Error: Required keys '%s' are missing from the object or 'sample'.", 
                 paste(missing_keys, collapse = ", ")))
  }
  
  # VALIDATION CHECK: Ensure no forbidden functions and disallowed variables coexist
  #----------------------------------------------------------------------------#
  
  # Forbidden functions and disallowed variables
  forbidden_functions <- c("n", "cur_group", "cur_group_id", "cur_group_rows", "cur_column")
  disallowed_variables <- c("abun", "rela", "norm")
  
  # Function to check for forbidden functions and disallowed variables in an expression
  check_expression_detailed <- function(expr, forbidden_functions, disallowed_variables) {
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
    check_result <- check_expression_detailed(expr_content, forbidden_functions, disallowed_variables)
    
    # If there is a violation, stop and show detailed error message
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
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  sample_info_data_merged <- object %>%
    purrr::map(\(x) {
      if(length(sample(x)) != 0){
        sample(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
    

  
  if (length(needed_abundance_keys) != 0){
    
    sample_abun_data_merged <- purrr::map(object, \(mgnet_obj){
      
      sample_abun_data <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                             taxa_id = taxa_id(mgnet_obj))
      
      for(abundance_key in needed_abundance_keys){
        sample_abun_data <- sample_abun_data %>%
          dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                             tibble::as_tibble(rownames = "sample_id") %>%
                             tidyr::pivot_longer(-sample_id, 
                                                 names_to = "taxa_id", 
                                                 values_to = abundance_key),
                           by = c("sample_id", "taxa_id"))
      }
      return(sample_abun_data)
    }) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      dplyr::left_join(sample_info_data_merged, 
                       by = c("mgnet","sample_id"))
  } 
  # END LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  # APPLY MUTATE
  #----------------------------------------------------------------------------#
  results <- sample_id(object) %>%
    purrr::map(\(x) tibble::tibble(sample_id = x)) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(quo_get_expr(expr))
    
    expr_name <- names(expressions)[i]
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      results <- sample_abun_data_merged %>%
        dplyr::group_by(!!!syms(.by)) %>%
        dplyr::mutate(!!(expr_name) := !!(expr)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("taxa_id", keys_abundance))) %>%
        dplyr::distinct() %>%
        dplyr::select(c("mgnet", "sample_id", tidyr::last_col())) %>%
        dplyr::right_join(results, by = c("mgnet","sample_id"))
        
      
    } else {
      
      # Process expressions not involving abundance-related variables
      results <- sample_info_data_merged %>%
        dplyr::group_by(!!!syms(.by)) %>%
        dplyr::mutate(!!(expr_name) := !!(expr)) %>%
        dplyr::ungroup() %>%
        dplyr::select(c("mgnet", "sample_id", tidyr::last_col())) %>%
        dplyr::right_join(results, by = c("mgnet","sample_id"))
      
    }
  }
  
  common_columns <- intersect(
    setdiff(names(sample_info_data_merged), c("sample_id", "mgnet")),
    setdiff(names(results), c("sample_id", "mgnet")))
  
  sample_info_data_merged <- sample_info_data_merged %>%
    dplyr::select(-dplyr::all_of(common_columns))
  
  sample_info_data_merged <- sample_info_data_merged %>%
    dplyr::left_join(results, dplyr::join_by("mgnet","sample_id"))
  
  sample_splitted <- sample_info_data_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  sample(object) <- sample_splitted
  return(object)
})