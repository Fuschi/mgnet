#' Filter Samples Data in mgnet or mgnetList Objects
#'
#' This function applies specified filtering criteria to samples data within an `mgnet` 
#' or `mgnetList` object. It intelligently uses reserved keywords to manage and manipulate
#' taxa data based on abundance metrics and samples metadata. The function is designed 
#' to integrate seamlessly with the tidyverse, supporting dplyr-style syntax and operations,
#' which facilitates complex data manipulation and filtering within a familiar framework.
#'
#' @param object An object of class \code{mgnet} or \code{mgnetList} containing the 
#'               data to be filtered.
#' @param ... Expressions that define the filtering criteria using variables available
#'            in the samples metadata and abundances of the object. These expressions can
#'            include tidyverse functions and operators, allowing for sophisticated 
#'            data manipulation directly within the function call.
#' @return Returns a modified \code{mgnet} or \code{mgnetList} object containing only 
#'         the samples that meet the specified filtering criteria. The method operates on each sample
#'         within an `mgnet` object individually. If the input is an `mgnetList`, the filtering 
#'         operation is applied independently to each `mgnet` object within the list, ensuring
#'         that results are consistent and isolated across different datasets or experimental conditions.
#'
#' @details
#' The function utilizes reserved keywords integral to the `mgnet` and `mgnetList` structures:
#'   - \code{abun}: raw abundance data.
#'   - \code{rela}: relative abundance data.
#'   - \code{norm}: normalized abundance data.
#'   - \code{sample_id}: A unique identifier for each sample.
#' Additionally, all other metadata columns within the sample data, in `meta` slot, can be used for filtering.
#' Internally, the method constructs a long-format data frame where each row represents a 
#' combination of sample and taxa data. This format facilitates complex filtering operations
#' across multiple dimensions of the data.
#'
#' @return Returns a modified \code{mgnet} or \code{mgnetList} object containing only 
#'         the samples that meet the specified filtering criteria.
#'
#' @export
#' @aliases filter_taxa_data,mgnet-method filter_taxa_data,mgnetList-method
#' @importFrom dplyr filter group_by ungroup arrange left_join select distinct
#' @importFrom rlang enquos eval_tidy
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect any_of
#' @importFrom purrr map map_lgl list_rbind
#' @importFrom methods slot
#' @importFrom magrittr %>%
setGeneric("filter_sample_data", function(object, ...) {standardGeneric("filter_sample_data")})

setMethod("filter_sample_data", "mgnet", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0 || ntaxa(object) == 0) {
    stop("Error: No samples or taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  # Check for missing non-abundance keys in sample
  if (length(needed_noabundance_keys) > 0) {
    missing_noabundance_keys <- setdiff(needed_noabundance_keys, sample_columns)
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
  long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                            taxa_id = taxa_id(object))
  
  for(abundance_key in needed_abundance_keys){
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        methods::slot(object, abundance_key) %>%
          tibble::as_tibble(rownames = "sample_id") %>%
          tidyr::pivot_longer(-sample_id, 
                              names_to = "taxa_id", 
                              values_to = abundance_key),
        by = c("sample_id", "taxa_id"))} 
  
  if(length(meta(object))!=0){
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(meta(object, .fmt = "tbl"), by = "sample_id")
  }
  
  # FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- long_mgnet %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_id, -tidyselect::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(sample_id, sample_id(object))) %>%
    dplyr::pull(sample_id)
    
  return(object[filtered_samples, ])
})

#------------------------------------------------------------------------------#
setMethod("filter_sample_data", "mgnetList", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0) || any(ntaxa(object) == 0)){
    stop("Error: No samples or taxa available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  
  
  # FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- sample_abun_data_merged %>%
    dplyr::group_by(mgnet, sample_id) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxa_id, -tidyselect::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::select(mgnet, sample_id)
  
  # SPLIT FILTERED
  #----------------------------------------------------------------------------#
  object <- purrr::imap(object, \(x, name){
    filtered_samples_name <- filtered_samples %>%
      dplyr::filter(mgnet == name) %>%
      dplyr::pull(sample_id)
    return(object[[name]][filtered_samples_name, ])
  })
  
  return(object)
})



#' Filter Taxa Data in mgnet or mgnetList Objects
#'
#' This function applies specified filtering criteria to taxa data within an `mgnet` 
#' or `mgnetList` object. It intelligently uses reserved keywords to manage and manipulate
#' taxa data based on abundance metrics and taxa metadata. The function is designed 
#' to integrate seamlessly with the tidyverse, supporting dplyr-style syntax and operations,
#' which facilitates complex data manipulation and filtering within a familiar framework.
#'
#' @param object An object of class \code{mgnet} or \code{mgnetList} containing the 
#'               data to be filtered.
#' @param ... Expressions that define the filtering criteria using variables available
#'            in the taxa metadata and abundances of the object. These expressions can
#'            include tidyverse functions and operators, allowing for sophisticated 
#'            data manipulation directly within the function call.
#' @return Returns a modified \code{mgnet} or \code{mgnetList} object containing only 
#'         the taxa that meet the specified filtering criteria. The method operates on each sample
#'         within an `mgnet` object individually. If the input is an `mgnetList`, the filtering 
#'         operation is applied independently to each `mgnet` object within the list, ensuring
#'         that results are consistent and isolated across different datasets or experimental conditions.
#'
#' @details
#' The function utilizes reserved keywords integral to the `mgnet` and `mgnetList` structures:
#'   - \code{abun}: raw abundance data.
#'   - \code{rela}: relative abundance data.
#'   - \code{norm}: normalized abundance data.
#'   - \code{taxa_id}: A unique identifier for each taxon.
#'   - \code{comm_id}: Community membership identifiers, used if the `comm` slot is present.
#' Additionally, all other metadata columns within the taxa data can be used for filtering.
#' Internally, the method constructs a long-format data frame where each row represents a 
#' combination of sample and taxa data. This format facilitates complex filtering operations
#' across multiple dimensions of the data.
#'
#' @return Returns a modified \code{mgnet} or \code{mgnetList} object containing only 
#'         the taxa that meet the specified filtering criteria.
#'
#' @export
#' @aliases filter_taxa_data,mgnet-method filter_taxa_data,mgnetList-method
#' @importFrom dplyr filter group_by ungroup arrange left_join select distinct
#' @importFrom rlang enquos eval_tidy
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect any_of
#' @importFrom purrr map map_lgl list_rbind
#' @importFrom methods slot
#' @importFrom magrittr %>%
setGeneric("filter_taxa_data", function(object, ...) {standardGeneric("filter_taxa_data")})

setMethod("filter_taxa_data", "mgnet", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0 || ntaxa(object) == 0) {
    stop("Error: No samples or taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
  long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                            taxa_id = taxa_id(object))
  
  for(abundance_key in needed_abundance_keys){
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        methods::slot(object, abundance_key) %>%
          tibble::as_tibble(rownames = "sample_id") %>%
          tidyr::pivot_longer(-sample_id, 
                              names_to = "taxa_id", 
                              values_to = abundance_key),
        by = c("sample_id", "taxa_id"))} 
  
  if(length(taxa(object))!=0){
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id")
  }
  
  # FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- long_mgnet %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sample_id, -tidyselect::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(taxa_id, taxa_id(object))) %>%
    dplyr::pull(taxa_id)
  
  return(object[, filtered_taxa])
})

#------------------------------------------------------------------------------#
setMethod("filter_taxa_data", "mgnetList", function(object, ...) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0) || any(ntaxa(object) == 0)){
    stop("Error: No samples or taxa available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
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
      dplyr::left_join(meta(object, .fmt = "tbl"), 
                       by = c("mgnet","sample_id"))
  }
  
  
  # FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- taxa_abun_data_merged %>%
    dplyr::group_by(mgnet, taxa_id) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sample_id, -tidyselect::any_of(keys_abundance)) %>%
    dplyr::distinct() %>%
    dplyr::select(mgnet, taxa_id)
  
  # SPLIT FILTERED
  #----------------------------------------------------------------------------#
  object <- purrr::imap(object, \(x, name){
    filtered_taxa_name <- filtered_taxa %>%
      dplyr::filter(mgnet == name) %>%
      dplyr::pull(taxa_id)
    return(object[[name]][, filtered_taxa_name])
  })
  
  return(object)
  
})