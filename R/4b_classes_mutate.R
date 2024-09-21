#' Modify and Augment `mgnet` Objects by Transforming the `sample` Slot
#'
#' This function dynamically manipulates the `sample` slot within `mgnet` or `mgnetList` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        The function targets the `sample` slot, which contains metadata for each sample.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the sample slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'sample_id' for `mgnet` objects and c('mgnet', 'sample_id') for `mgnetList` objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by `.by`. Usage of 'taxa_id' as a grouping variable is strictly prohibited to maintain
#'        data consistency and avoid misinterpretation.
#'
#' @details The function is designed to integrate seamlessly with the `tidyverse`, allowing users
#'          to utilize familiar and potent data manipulation verbs such as `mutate`, `filter`.
#'          It supports using any `tidyverse`-compatible expressions, including conditional operations,
#'          summarizations, and complex transformations involving both abundance and metadata fields.
#'          This flexibility makes it particularly useful for ecological and biological data analysis,
#'          where combining different data types and conditions is common.
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **sample_id**: An essential identifier used to uniquely reference individual samples within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnetList` object with updated `sample` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_sample,mgnet-method mutate_sample,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate
#' @importFrom tidyr expand_grid any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_sample", function(object, ..., .by) {standardGeneric("mutate_sample")})

setMethod("mutate_sample", "mgnet", function(object, ..., .by = "sample_id") {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "sample_id"
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
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
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  if(length(sample(object)) != 0){
    sample_info_mutated <- sample(object, .fmt = "tbl")
  } else {
    sample_info_mutated <- tibble::tibble(sample_id = sample_id(object))
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])

    if(any(expr_vars %in% c("abun","rela","norm"))){
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      sample_abun_data <- tidyr::expand_grid(sample_id = sample_id(object),
                                             taxa_id = taxa_id(object))
      
      for(abundance_key in intersect(expr_vars, c("abun","rela", "norm"))){
        sample_abun_data <- sample_abun_data %>%
          dplyr::left_join(
            methods::slot(object, abundance_key) %>%
              tibble::as_tibble(rownames = "sample_id") %>%
              tidyr::pivot_longer(-sample_id, 
                                  names_to = "taxa_id", 
                                  values_to = abundance_key),
            by = c("sample_id", "taxa_id"))} 
      
      sample_abun_data <- sample_abun_data %>%
        left_join(sample_info_mutated, by = "sample_id")
      
      sample_info_mutated <- sample_abun_data %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::arrange(match(sample_id, sample_id(object))) 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      sample_info_mutated <- sample_info_mutated %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(match(sample_id, sample_id(object))) 
    }
  } # loop expressions
  
  sample(object) <- tibble::column_to_rownames(sample_info_mutated, "sample_id")
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
  expressions <- rlang::enquos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","sample_id")
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
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

  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_sample_mutated_merged <- object %>%
    purrr::map(\(x) {
      if(length(sample(x)) != 0){
        sample(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      sample_abun_data_merged <- purrr::map(object, \(mgnet_obj){
        
        sample_abun_data <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                               taxa_id = taxa_id(mgnet_obj))
        
        for(abundance_key in intersect(expr_vars, c("abun","rela", "norm"))){
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
        dplyr::left_join(info_sample_mutated_merged, 
                         by = c("mgnet","sample_id"))
      
      info_sample_mutated_merged <- sample_abun_data_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      info_sample_mutated_merged <- info_sample_mutated_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  } # loop expressions
  
  info_sample_mutated_splitted <- info_sample_mutated_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  sample(object) <- info_sample_mutated_splitted
  return(object)
})


#' Modify and Augment `mgnet` Objects by Transforming the `taxa` Slot
#'
#' This function dynamically manipulates the `taxa` slot within `mgnet` or `mgnetList` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        The function targets the `taxa` slot, which contains metadata for each taxa
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the sample slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'taxa_id' for `mgnet` objects and c('mgnet', 'taxa_id') for `mgnetList` objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by `.by`. Usage of 'sample_id' as a grouping variable is strictly prohibited to maintain
#'        data consistency and avoid misinterpretation.
#'
#' @details The function is designed to integrate seamlessly with the `tidyverse`, allowing users
#'          to utilize familiar and potent data manipulation verbs such as `mutate`, `filter`.
#'          It supports using any `tidyverse`-compatible expressions, including conditional operations,
#'          summarizations, and complex transformations involving both abundance and metadata fields.
#'          This flexibility makes it particularly useful for ecological and biological data analysis,
#'          where combining different data types and conditions is common.
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **taxa_id**: An essential identifier used to uniquely reference individual taxon within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnetList` object with updated `taxa` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_taxa,mgnet-method mutate_taxa,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate
#' @importFrom tidyr expand_grid any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_taxa", function(object, ..., .by) {standardGeneric("mutate_taxa")})

setMethod("mutate_taxa", "mgnet", function(object, ..., .by = "taxa_id") {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (ntaxa(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "taxa_id"
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
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
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  taxa_info_mutated <- tibble::tibble(taxa_id = taxa_id(object))
  if(length(taxa(object))!=0){
    taxa_info_mutated <- taxa_info_mutated %>% 
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id")
  }
  if(length(comm(x))!=0){
    taxa_info_mutated <- taxa_info_mutated %>% 
      dplyr::left_join(comm_id(object, .fmt = "tbl"), by = "taxa_id")
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      taxa_abun_data <- tidyr::expand_grid(sample_id = sample_id(object),
                                             taxa_id = taxa_id(object))
      
      for(abundance_key in intersect(expr_vars, c("abun","rela", "norm"))){
        taxa_abun_data <- taxa_abun_data %>%
          dplyr::left_join(
            methods::slot(object, abundance_key) %>%
              tibble::as_tibble(rownames = "sample_id") %>%
              tidyr::pivot_longer(-sample_id, 
                                  names_to = "taxa_id", 
                                  values_to = abundance_key),
            by = c("sample_id", "taxa_id"))} 
      
      taxa_abun_data <- taxa_abun_data %>%
        left_join(taxa_info_mutated, by = "taxa_id")
      
      taxa_info_mutated <- taxa_abun_data %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::arrange(match(taxa_id, taxa_id(object))) 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      taxa_info_mutated <- taxa_info_mutated %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(match(taxa_id, taxa_id(object))) 
    }
  } # loop expressions
  
  taxa(object) <- tibble::column_to_rownames(taxa_info_mutated, "taxa_id")
  return(object)
  
})


setMethod("mutate_taxa", "mgnetList", function(object, ..., .by = c("mgnet", "taxa_id")) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Ensure there are samples to process
  if (any(ntaxa(object) == 0)) {
    stop("Error: No samples available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","taxa_id")
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
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
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_taxa_mutated_merged <- object %>%
    purrr::map(\(x) {
      y <- tibble::tibble(taxa_id = taxa_id(x))
      if(length(taxa(object))!=0){
        y <- y %>% dplyr::left_join(taxa(x, .fmt = "tbl"), by = "taxa_id")
      }
      if(length(comm(x))!=0){
        y <- y %>% dplyr::left_join(comm_id(x, .fmt = "tbl"), by = "taxa_id")
      }
      return(y)
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      taxa_abun_data_merged <- purrr::map(object, \(mgnet_obj){
        
        taxa_abun_data <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                             taxa_id = taxa_id(mgnet_obj))
        
        for(abundance_key in intersect(expr_vars, c("abun","rela", "norm"))){
          taxa_abun_data <- taxa_abun_data %>%
            dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                               tibble::as_tibble(rownames = "sample_id") %>%
                               tidyr::pivot_longer(-sample_id, 
                                                   names_to = "taxa_id", 
                                                   values_to = abundance_key),
                             by = c("sample_id", "taxa_id"))
        }
        return(taxa_abun_data)
      }) %>%
        purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
        purrr::list_rbind() %>%
        dplyr::left_join(info_taxa_mutated_merged, 
                         by = c("mgnet","taxa_id"))
      
      info_taxa_mutated_merged <- taxa_abun_data_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      info_taxa_mutated_merged <- info_taxa_mutated_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  } # loop expressions
  
  info_taxa_mutated_splitted <- info_taxa_mutated_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  taxa(object) <- info_taxa_mutated_splitted
  return(object)
})