#------------------------------------------------------------------------------#
#' Modify and Augment Sample Data in `mgnet` and `mgnetList` Objects
#'
#' This function dynamically manipulates both the sample metadata and abundance data within 
#' `mgnet` and `mgnetList` objects, applying user-defined transformations. It leverages the `tidyverse` tools 
#' such as `dplyr` to enable flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object. For `mgnet`, the function modifies sample 
#'        metadata and abundance data within the object. For `mgnetList`, the same transformation 
#'        is applied to each `mgnet` object in the list.
#' @param ... Dynamic expressions or functions to be applied to the sample data and metadata.
#'        These expressions can include both abundance-related fields (e.g., 'abun', 'rela', 'norm') 
#'        and metadata fields.
#'
#' @section Reserved Keywords:
#' Reserved keywords such as `sample_id`, `taxa_id`, `comm_id`, `abun`, `rela`, `norm`, and `mgnet`
#' are essential to the structure of the `mgnet` and `mgnetList` objects. They should not be 
#' overwritten or modified by user-defined transformations. Any attempt to modify these keywords 
#' will result in an error.
#'
#' @return Returns the modified `mgnet` or `mgnetList` object with the transformed `meta` slot 
#'         reflecting the applied transformations. All other structures within the object remain 
#'         unchanged.
#'
#' @export
#' @aliases mutate_sample_data,mgnet-method mutate_sample_data,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup arrange left_join group_by select distinct
#' @importFrom rlang enquos syms eval_tidy quo_name quo_get_expr
#' @importFrom tibble column_to_rownames tibble add_column as_tibble
#' @importFrom purrr list_rbind map imap map_lgl
#' @importFrom methods slot
#' @importFrom tidyr any_of pivot_longer expand_grid
setGeneric("mutate_sample_data", function(object, ...) {standardGeneric("mutate_sample_data")})

setMethod("mutate_sample_data", "mgnet", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0 || ntaxa(object) == 0) {
    stop("Error: No samples or taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")
  
  # Get sample metadata variables
  sample_columns <- meta_vars(object)
  
  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Check if at least one abundance key is present
  if (length(needed_abundance_keys) == 0) {
    stop("Error: At least one abundance-related key ('abun', 'rela', 'norm') must be present.")
  }
  
  # Check each abundance key for presence
  missing_abundances <- c()
  for (key in keys_abundance) {
    if (key %in% keys_required && length(methods::slot(object, key)) == 0) {
      missing_abundances <- c(missing_abundances, key)
    }
  }
  
  # Check for missing non-abundance keys in sample
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundances <- setdiff(needed_noabundance_keys, sample_columns)
  }
  
  # Generate error message if any keys are missing
  if (length(missing_abundances) > 0 || length(missing_noabundances) > 0) {
    error_message <- "Error: The following required variables in the expressions are missing:\n"
    
    # Add missing abundance-related keys
    if (length(missing_abundances) > 0) {
      error_message <- paste0(error_message, "- Abundance-related: ", 
                              paste(missing_abundances, collapse = ", "), "\n")
    }
    
    # Add missing sample metadata keys
    if (length(missing_noabundances) > 0) {
      error_message <- paste0(error_message, "- Sample metadata: ", 
                              paste(missing_noabundances, collapse = ", "), "\n")
    }
    
    stop(error_message)
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT DATA
  #----------------------------------------------------------------------------#
  sample_mutated_info <- tidyr::expand_grid(sample_id = sample_id(object),
                                            taxa_id = taxa_id(object))
  
  for(abundance_key in needed_abundance_keys){
    sample_mutated_info <- sample_mutated_info %>%
      dplyr::left_join(
        methods::slot(object, abundance_key) %>%
          tibble::as_tibble(rownames = "sample_id") %>%
          tidyr::pivot_longer(-sample_id, 
                              names_to = "taxa_id", 
                              values_to = abundance_key),
        by = c("sample_id", "taxa_id"))} 
  
  if(length(meta(object))!=0){
    sample_mutated_info <- sample_mutated_info %>%
      dplyr::left_join(meta(object, .fmt = "tbl"), by = "sample_id")
  }
  
  # MUTATE
  #----------------------------------------------------------------------------#
  sample_mutated_info <- sample_mutated_info %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_id, -tidyr::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(sample_id, sample_id(object))) %>%
    tibble::column_to_rownames("sample_id")
    
  meta(object) <- sample_mutated_info
  return(object)
})

#------------------------------------------------------------------------------#
setMethod("mutate_sample_data", "mgnetList", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0) || any(ntaxa(object) == 0)){
    stop("Error: No samples or taxa available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")
  
  # Get sample metadata variables
  sample_columns <- unique(unlist(meta_vars(object)))
  
  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Prepare to check for missing keys across all mgnet objects in mgnetList
  missing_abundances <- c()
  missing_noabundances <- c()
  
  # Check if at least one abundance key is present in this mgnet object
  if (length(intersect(needed_abundance_keys, keys_abundance)) == 0) {
    stop(sprintf("Error: At least one abundance-related key ('abun', 'rela', 'norm') must be present in mgnet object %d.", i))
  }
  
  # Iterate over each mgnet object in the list
  for (i in seq_along(object)) {
    mgnet_obj <- object[[i]]
    
    # Check each abundance key for presence in the current mgnet object
    for (key in keys_abundance) {
      if (key %in% keys_required && length(methods::slot(mgnet_obj, key)) == 0) {
        missing_abundances <- c(missing_abundances, sprintf("%s in mgnet object %d", key, i))
      }
    }
  }
  
  # Generate error message if any specific keys are missing
  if (length(missing_abundances) > 0 || length(missing_noabundances) > 0) {
    error_message <- "Error: The following required keys are missing:\n"
    
    # Add missing abundance-related keys
    if (length(missing_abundances) > 0) {
      error_message <- paste0(error_message, "- Abundance-related: ", 
                              paste(missing_abundances, collapse = ", "), "\n")
    }
    
    # Add missing sample metadata keys
    if (length(missing_noabundances) > 0) {
      error_message <- paste0(error_message, "- Sample metadata: ", 
                              paste(missing_noabundances, collapse = ", "), "\n")
    }
    
    stop(error_message)
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT DATA
  #----------------------------------------------------------------------------#
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
    purrr::list_rbind() 
  
  if(length(needed_noabundance_keys) != 0){
    sample_abun_data_merged <- sample_abun_data_merged %>%
      dplyr::left_join(meta(object, .fmt = "tbl"), 
                       by = c("mgnet","sample_id"))
    }
  
  
  # MUTATE
  #----------------------------------------------------------------------------#
  sample_abun_data_merged <- sample_abun_data_merged %>%
    dplyr::group_by(mgnet, sample_id) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_id, -tidyr::any_of(keys_abundance)) %>%
    dplyr::distinct() 
  
  # SPLIT MUTATED
  #----------------------------------------------------------------------------#
  sample_abun_data_splitted <- sample_abun_data_merged %>%
  base::split(.[, "mgnet"]) %>%
  purrr::imap(\(x,y){
    dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
  }) %>%
  purrr::map(\(x){
    x %>% dplyr::select(-"mgnet") %>%
      tibble::column_to_rownames("sample_id")
  })
  
  meta(object) <- sample_abun_data_splitted
  return(object)
})


#------------------------------------------------------------------------------#
#' Modify and Augment Taxa Data in `mgnet` and `mgnetList` Objects
#'
#' This function dynamically manipulates both the taxa metadata and abundance data within 
#' `mgnet` and `mgnetList` objects, applying user-defined transformations. It leverages the `tidyverse` tools 
#' such as `dplyr` to enable flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object. For `mgnet`, the function modifies taxa 
#'        metadata and abundance data within the object. For `mgnetList`, the same transformation 
#'        is applied to each `mgnet` object in the list.
#' @param ... Dynamic expressions or functions to be applied to the taxa data and metadata.
#'        These expressions can include both abundance-related fields (e.g., 'abun', 'rela', 'norm') 
#'        and metadata fields.
#'
#' @section Reserved Keywords:
#' Reserved keywords such as `sample_id`, `taxa_id`, `comm_id`, `abun`, `rela`, `norm`, and `mgnet`
#' are essential to the structure of the `mgnet` and `mgnetList` objects. They should not be 
#' overwritten or modified by user-defined transformations. Any attempt to modify these keywords 
#' will result in an error.
#'
#' @return Returns the modified `mgnet` or `mgnetList` object with the transformed `meta` slot 
#'         reflecting the applied transformations. All other structures within the object remain 
#'         unchanged.
#'
#' @export
#' @aliases mutate_taxa_data,mgnet-method mutate_taxa_data,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup arrange left_join group_by select distinct
#' @importFrom rlang enquos syms eval_tidy quo_name quo_get_expr
#' @importFrom tibble column_to_rownames tibble add_column as_tibble
#' @importFrom purrr list_rbind map imap map_lgl
#' @importFrom methods slot
#' @importFrom tidyr any_of pivot_longer expand_grid
setGeneric("mutate_taxa_data", function(object, ...) {standardGeneric("mutate_taxa_data")})

setMethod("mutate_taxa_data", "mgnet", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0 || ntaxa(object) == 0) {
    stop("Error: No samples or taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")
  
  # Get taxa metadata variables
  taxa_columns <- taxa_vars(object)
  
  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Prepare to check for missing keys
  missing_abundances <- c()
  missing_noabundances <- c()
  
  # Check if at least one abundance key is present
  if (length(needed_abundance_keys) == 0) {
    stop("Error: At least one abundance-related key ('abun', 'rela', 'norm') must be present.")
  }
  
  # Check each abundance key for presence
  for (key in keys_abundance) {
    if (key %in% keys_required && length(methods::slot(object, key)) == 0) {
      missing_abundances <- c(missing_abundances, key)
    }
  }
  
  # Check for missing non-abundance keys in taxa
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundance_keys <- setdiff(needed_noabundance_keys, taxa_columns)
    if (length(missing_noabundance_keys) > 0) {
      missing_noabundances <- c(missing_noabundances, missing_noabundance_keys)
    }
  }
  
  # Generate error message if any keys are missing
  if (length(missing_abundances) > 0 || length(missing_noabundances) > 0) {
    error_message <- "Error: The following required variables in the expressions are missing:\n"
    
    # Add missing abundance-related keys
    if (length(missing_abundances) > 0) {
      error_message <- paste0(error_message, "- Abundance-related: ", 
                              paste(missing_abundances, collapse = ", "), "\n")
    }
    
    # Add missing taxa metadata keys
    if (length(missing_noabundances) > 0) {
      error_message <- paste0(error_message, "- taxa metadata: ", 
                              paste(missing_noabundances, collapse = ", "), "\n")
    }
    
    stop(error_message)
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT DATA
  #----------------------------------------------------------------------------#
  taxa_mutated_info <- tidyr::expand_grid(sample_id = sample_id(object),
                                            taxa_id = taxa_id(object))
  
  for(abundance_key in needed_abundance_keys){
    taxa_mutated_info <- taxa_mutated_info %>%
      dplyr::left_join(
        methods::slot(object, abundance_key) %>%
          tibble::as_tibble(rownames = "sample_id") %>%
          tidyr::pivot_longer(-sample_id, 
                              names_to = "taxa_id", 
                              values_to = abundance_key),
        by = c("sample_id", "taxa_id"))} 
  
  if(length(taxa(object))!=0){
    taxa_mutated_info <- taxa_mutated_info %>%
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id")
  }
  
  # MUTATE
  #----------------------------------------------------------------------------#
  taxa_mutated_info <- taxa_mutated_info %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_id, -tidyr::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(taxa_id, taxa_id(object))) %>%
    tibble::column_to_rownames("taxa_id")
  
  taxa(object) <- taxa_mutated_info
  return(object)
})

#------------------------------------------------------------------------------#
setMethod("mutate_taxa_data", "mgnetList", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0) || any(ntaxa(object) == 0)){
    stop("Error: No samples or taxa available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Define abundance-related keys
  keys_abundance <- c("abun", "rela", "norm")
  
  # Get taxa metadata variables
  taxa_columns <- unique(unlist(meta_vars(object)))
  
  # Determine which abundance keys are needed
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Prepare to check for missing keys across all mgnet objects in mgnetList
  missing_abundances <- c()
  missing_noabundances <- c()
  
  # Check if at least one abundance key is present in this mgnet object
  if (length(intersect(needed_abundance_keys, keys_abundance)) == 0) {
    stop(sprintf("Error: At least one abundance-related key ('abun', 'rela', 'norm') must be present in mgnet object %d.", i))
  }
  
  # Iterate over each mgnet object in the list
  for (i in seq_along(object)) {
    mgnet_obj <- object[[i]]
    
    # Check each abundance key for presence in the current mgnet object
    for (key in keys_abundance) {
      if (key %in% keys_required && length(methods::slot(mgnet_obj, key)) == 0) {
        missing_abundances <- c(missing_abundances, sprintf("%s in mgnet object %d", key, i))
      }
    }
  }
  
  # Generate error message if any specific keys are missing
  if (length(missing_abundances) > 0 || length(missing_noabundances) > 0) {
    error_message <- "Error: The following required keys are missing:\n"
    
    # Add missing abundance-related keys
    if (length(missing_abundances) > 0) {
      error_message <- paste0(error_message, "- Abundance-related: ", 
                              paste(missing_abundances, collapse = ", "), "\n")
    }
    
    # Add missing taxa metadata keys
    if (length(missing_noabundances) > 0) {
      error_message <- paste0(error_message, "- taxa metadata: ", 
                              paste(missing_noabundances, collapse = ", "), "\n")
    }
    
    stop(error_message)
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT DATA
  #----------------------------------------------------------------------------#
  taxa_abun_data_merged <- purrr::map(object, \(mgnet_obj){
    
    taxa_abun_data <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                           taxa_id = taxa_id(mgnet_obj))
    
    for(abundance_key in needed_abundance_keys){
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
    purrr::list_rbind() 
  
  if(length(needed_noabundance_keys) != 0){
    taxa_abun_data_merged <- taxa_abun_data_merged %>%
      dplyr::left_join(taxa(object, .fmt = "tbl"), 
                       by = c("mgnet","taxa_id"))
  }
  
  
  # MUTATE
  #----------------------------------------------------------------------------#
  taxa_abun_data_merged <- taxa_abun_data_merged %>%
    dplyr::group_by(mgnet, taxa_id) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sample_id, -tidyr::any_of(keys_abundance)) %>%
    dplyr::distinct() 
  
  # SPLIT MUTATED
  #----------------------------------------------------------------------------#
  taxa_abun_data_splitted <- taxa_abun_data_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  taxa(object) <- taxa_abun_data_splitted
  return(object)
})