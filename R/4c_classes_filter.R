#' Filter `mgnet` or `mgnetList` Objects Based on Sample and Abundance Data
#'
#' This function filters entire `mgnet` or `mgnetList` objects, including their sample and abundance data,
#' based on user-specified conditions. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        This function operates on the entire structure of these objects, allowing for the filtering
#'        of samples based on both metadata in the sample field and various forms of abundance data.
#' @param ... Dynamic expressions or functions to apply for filtering.
#'        Users can specify conditions involving both metadata and abundance data, facilitating
#'        complex and targeted filtering strategies that are essential in data-driven studies.
#' @param .by Optional; a character vector specifying the columns to group data by before filtering.
#'        Defaults to 'sample_id' for `mgnet` objects and c('mgnet', 'sample_id') for `mgnetList` objects.
#'        Proper grouping is essential for ensuring that filters are applied contextually and correctly.
#'        The use of 'taxa_id' as a grouping variable is prohibited to maintain clarity and data integrity.
#'
#' @details Leveraging the `tidyverse` functionalities, this function allows for the integration of
#'          sophisticated data manipulation techniques. It supports conditional operations, group-based
#'          filtering, and data summarization, making it highly effective for analyses where detailed
#'          sample selection and data refinement are required. The ability to filter using both sample
#'          metadata and abundance data offers unparalleled flexibility and power in data handling.
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **sample_id**: An essential identifier used to uniquely reference individual samples within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list supporting sophisticated multi-dataset management.
#'
#' @return Returns the `mgnet` or `mgnetList` object with updated content reflecting the applied filters.
#'         This includes any transformations to sample metadata and adjustments based on abundance data.
#'         The structure outside the targeted data fields remains unchanged, preserving the integrity of the objects.
#'
#' @export
#' @aliases filter_sample,mgnet-method filter_sample,mgnetList-method
#' @importFrom dplyr filter group_by ungroup distinct relocate semi_join
#' @importFrom tidyr expand_grid any_of pivot_longer
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr reduce map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("filter_sample", function(object, ..., .by) {standardGeneric("filter_sample")})

setMethod("filter_sample", "mgnet", function(object, ..., .by = "sample_id") {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Ensure there are samples to process
  if (nsample(object) == 0) {
    stop("Error: No samples available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Check the variables needed
  validate_required_variables(object, keys_required, "sample")
  
  # Store needed abundances keys
  keys_abundance <- c("abun", "rela", "norm")
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "sample_id"
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  if(length(meta(object)) == 0){
    sample_info_data = tibble::tibble(sample_id = sample_id(object))
  } else {
    sample_info_data = meta(object, .fmt = "tbl")
  }
  
  if(length(needed_abundance_keys) != 0){
    
    long_abun <- tidyr::expand_grid(sample_id = sample_id(object),
                                           taxa_id = taxa_id(object))
    
    for(abundance_key in needed_abundance_keys){
      long_abun <- long_abun %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id, 
                                names_to = "taxa_id", 
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
    
    if(length(needed_noabundance_keys) != 0){
      long_abun <- long_abun %>%
        dplyr::left_join(meta(object, .fmt = "tbl"), by = "sample_id")
    }
    
  }
  # END LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- list()
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(rlang::quo_get_expr(expr))
    
    #expr_name <- names(expressions)[i]
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      filtered_samples[[i]] <- long_abun %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("taxa_id", keys_abundance))) %>%
        dplyr::distinct() %>%
        dplyr::pull("sample_id")

    } else {
      
      # Process expressions not involving abundance-related variables
      filtered_samples[[i]] <- sample_info_data %>%
          dplyr::group_by(!!!rlang::syms(.by)) %>%
          dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
          dplyr::ungroup() %>%
          dplyr::pull("sample_id")

    }
  }
    
  filtered_samples <- purrr::reduce(filtered_samples, intersect) 
  filtered_samples <- filtered_samples[order(which(sample_id(object) %in% filtered_samples))]

  return(object[filtered_samples, ])
})


#------------------------------------------------------------------------------#
setMethod("filter_sample", "mgnetList", function(object, ..., .by = c("mgnet", "sample_id")) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0)) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Check the variables needed
  lapply(object, \(x){
    validate_required_variables(x, keys_required, "sample")})
  
  # Store needed abundances keys
  keys_abundance <- c("abun","rela","norm")
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","sample_id")
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  sample_info_data_merged <- object %>%
    purrr::map(\(x) {
      if(length(meta(x)) != 0){
        meta(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  if (length(needed_abundance_keys) != 0){
    
    long_abun_merged <- purrr::map(object, \(mgnet_obj){
      
      long_abun <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                             taxa_id = taxa_id(mgnet_obj))
      
      for(abundance_key in needed_abundance_keys){
        long_abun <- long_abun %>%
          dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                             tibble::as_tibble(rownames = "sample_id") %>%
                             tidyr::pivot_longer(-sample_id, 
                                                 names_to = "taxa_id", 
                                                 values_to = abundance_key),
                           by = c("sample_id", "taxa_id"))
      }
      return(long_abun)
    }) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      dplyr::left_join(sample_info_data_merged, 
                       by = c("mgnet","sample_id"))
  } 
  # END LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- list()
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(quo_get_expr(expr))
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      filtered_samples[[i]] <- long_abun_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyr::all_of(c("mgnet", "sample_id")))
      
    } else {
      
      # Process expressions not involving abundance-related variables
      filtered_samples[[i]] <- sample_info_data_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  }
  
  filtered_samples <- purrr::reduce(filtered_samples, dplyr::semi_join, by = c("mgnet", "sample_id")) 
  
  for(i in names(object)){
    filtered_i <- dplyr::filter(filtered_samples, mgnet == i) %>% dplyr::pull("sample_id")
    filtered_i <- filtered_i[order(which(sample_id(object[[i]]) %in% filtered_i))]
    object[[i]] <- object[[i]][filtered_i, ]
  }
  
  validObject(object)
  return(object)
})




#' Filter `mgnet` or `mgnetList` Objects Based on Taxa and Abundance Data
#'
#' This function filters entire `mgnet` or `mgnetList` objects, including their taxa and abundance data,
#' based on user-specified conditions. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        This function operates on the entire structure of these objects, allowing for the filtering
#'        of taxa based on both metadata in the taxa field and various forms of abundance data.
#' @param ... Dynamic expressions or functions to apply for filtering.
#'        Users can specify conditions involving both metadata and abundance data, facilitating
#'        complex and targeted filtering strategies that are essential in data-driven studies.
#' @param .by Optional; a character vector specifying the columns to group data by before filtering.
#'        Defaults to 'taxa_id' for `mgnet` objects and c('mgnet', 'taxa_id') for `mgnetList` objects.
#'        Proper grouping is essential for ensuring that filters are applied contextually and correctly.
#'        The use of 'sample_id' as a grouping variable is prohibited to maintain clarity and data integrity.
#'
#' @details Leveraging the `tidyverse` functionalities, this function allows for the integration of
#'          sophisticated data manipulation techniques. It supports conditional operations, group-based
#'          filtering, making it highly effective for analyses where detailed
#'          taxa selection and data refinement are required. 
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **taxa_id**: An essential identifier used to uniquely reference individual taxa within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list supporting sophisticated multi-dataset management.
#'
#' @return Returns the `mgnet` or `mgnetList` object with updated content reflecting the applied filters.
#'         This includes any transformations to taxa metadata and adjustments based on abundance data.
#'         The structure outside the targeted data fields remains unchanged, preserving the integrity of the objects.
#'
#' @export
#' @aliases filter_taxa,mgnet-method filter_taxa,mgnetList-method
#' @importFrom dplyr filter group_by ungroup distinct relocate semi_join
#' @importFrom tidyr expand_grid any_of pivot_longer
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr reduce map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("filter_taxa", function(object, ..., .by) {standardGeneric("filter_taxa")})

setMethod("filter_taxa", "mgnet", function(object, ..., .by = "taxa_id") {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (ntaxa(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Check the variables needed
  validate_required_variables(object, keys_required, "taxa")
  
  # Store needed abundances keys
  keys_abundance <- c("abun","rela","norm")
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "taxa_id"
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE LONG FORMAT taxa DATA
  #----------------------------------------------------------------------------#
  
  if(length(taxa(object)) == 0){
    taxa_info_data = tibble::tibble(taxa_id = taxa_id(object))
  } else {
    taxa_info_data = taxa(object, .fmt = "tbl")
  }
  
  taxa_info_data <- tibble::tibble(taxa_id = taxa_id(object))
  if(length(taxa(object))!=0){
    taxa_info_data <- taxa_info_data %>% 
      dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id")
  }
  if(length(comm(object))!=0){
    taxa_info_data <- taxa_info_data %>% 
      dplyr::left_join(comm_id(object, .fmt = "tbl"), by = "taxa_id")
  }
  
  
  if(length(needed_abundance_keys) != 0){
    
    taxa_abun_data <- tidyr::expand_grid(sample_id = sample_id(object),
                                           taxa_id = taxa_id(object))
    
    for(abundance_key in needed_abundance_keys){
      taxa_abun_data <- taxa_abun_data %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id, 
                                names_to = "taxa_id", 
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
    
    if(length(taxa(object)) != 0){
      taxa_abun_data <- taxa_abun_data %>%
        dplyr::left_join(taxa(object, .fmt = "tbl"),
                         by = "taxa_id")
    }
    
  }
  # END LONG FORMAT SAMPLE DATA
  #----------------------------------------------------------------------------#
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- list()
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(rlang::quo_get_expr(expr))
    
    #expr_name <- names(expressions)[i]
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      filtered_taxa[[i]] <- taxa_abun_data %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("sample_id", keys_abundance))) %>%
        dplyr::distinct() %>%
        dplyr::pull("taxa_id")
      
    } else {
      
      # Process expressions not involving abundance-related variables
      filtered_taxa[[i]] <- taxa_info_data %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::pull("taxa_id")
      
    }
  }
  
  filtered_taxa <- purrr::reduce(filtered_taxa, intersect) 
  filtered_taxa <- filtered_taxa[order(which(taxa_id(object) %in% filtered_taxa))]
  
  return(object[, filtered_taxa])
})


setMethod("filter_taxa", "mgnetList", function(object, ..., .by = c("mgnet", "taxa_id")) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(ntaxa(object) == 0)) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Check the variables needed
  validate_required_variables(object, keys_required, "taxa")
  
  # Store needed abundances keys
  keys_abundance <- c("abun","rela","norm")
  needed_abundance_keys <- intersect(keys_required, keys_abundance)
  needed_noabundance_keys <- setdiff(keys_required, keys_abundance)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","taxa_id")
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # CREATE LONG FORMAT taxa DATA
  #----------------------------------------------------------------------------#

  taxa_info_data_merged <- object %>%
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
  
  if (length(needed_abundance_keys) != 0){
    
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
      purrr::list_rbind() %>%
      dplyr::left_join(taxa_info_data_merged, 
                       by = c("mgnet","taxa_id"))
  } 
  # END LONG FORMAT taxa DATA
  #----------------------------------------------------------------------------#
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- list()
  
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    expr_vars <- all.vars(quo_get_expr(expr))
    
    if (any(expr_vars %in% keys_abundance)) {
      
      # Process expressions involving abundance-related variables
      filtered_taxa[[i]] <- taxa_abun_data_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyr::all_of(c("mgnet", "taxa_id")))
      
    } else {
      
      # Process expressions not involving abundance-related variables
      filtered_taxa[[i]] <- taxa_info_data_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  }
  
  filtered_taxa <- purrr::reduce(filtered_taxa, dplyr::semi_join, by = c("mgnet", "taxa_id")) 
  
  for(i in names(object)){
    filtered_i <- dplyr::filter(filtered_taxa, mgnet == i) %>% dplyr::pull("taxa_id")
    filtered_i <- filtered_i[order(which(taxa_id(object[[i]]) %in% filtered_i))]
    object[[i]] <- object[[i]][, filtered_i]
  }
  
  validObject(object)
  return(object)
})