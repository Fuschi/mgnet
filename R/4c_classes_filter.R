#' Filter `mgnet` or `mgnetList` Objects Based on Sample Metadata and Abundances Criteria
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
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'sample_id' for mgnet objects and c('mgnet', 'sample_id') for mgnetList objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by .by. If you do not wish to group the data, set .by to NULL.
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
#' @aliases filter_meta,mgnet-method filter_meta,mgnetList-method
#' @importFrom dplyr filter group_by ungroup distinct relocate semi_join
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom tidyselect any_of all_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr reduce map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("filter_meta", function(object, ..., .by = NULL) {standardGeneric("filter_meta")})

setMethod("filter_meta", "mgnet", function(object, ..., .by = NULL) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are sample to process
  if(miss_sample(object)) stop("Error: No sample available in the 'mgnet' object.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- "sample_id"
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#

  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  metadata <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- apply_filter_verb(metadata, long_abun, expressions, .by, "sample")
  
  filtered_samples <- lapply(filtered_samples, \(x) dplyr::pull(x, sample_id))
  filtered_samples <- purrr::reduce(filtered_samples, intersect) 
  filtered_samples <- filtered_samples[order(which(sample_id(object) %in% filtered_samples))]

  return(object[filtered_samples, ])
})


#------------------------------------------------------------------------------#
setMethod("filter_meta", "mgnetList", function(object, ..., .by = c("mgnet", "sample_id")) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are sample to process
  if(miss_sample(object, "any")) stop("Error: No sample available in at least one of the mgnet objects.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- c("mgnet", "sample_id")
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  metadata <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_samples <- apply_filter_verb(metadata, long_abun, expressions, .by, "sample") 
  
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
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'taxa_id' for mgnet objects and c('mgnet', 'taxa_id') for mgnetList objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by .by. If you do not wish to group the data, set .by to NULL.
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
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom tidyselect any_of all_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr reduce map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("filter_taxa", function(object, ..., .by = NULL) {standardGeneric("filter_taxa")})

setMethod("filter_taxa", "mgnet", function(object, ..., .by = NULL) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are taxa to process
  if(miss_taxa(object)) stop("Error: No taxa available in the 'mgnet' object.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- "taxa_id"
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  metadata <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- apply_filter_verb(metadata, long_abun, expressions, .by, "taxa")
  
  filtered_taxa <- lapply(filtered_taxa, \(x) dplyr::pull(x, taxa_id))
  filtered_taxa <- purrr::reduce(filtered_taxa, intersect) 
  filtered_taxa <- filtered_taxa[order(which(taxa_id(object) %in% filtered_taxa))]
  return(object[, filtered_taxa])
})


setMethod("filter_taxa", "mgnetList", function(object, ..., .by = NULL) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are taxa to process
  if(miss_taxa(object, "any")) stop("Error: No taxa available in at least one of the mgnet objects.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- c("mgnet", "taxa_id")

  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  metadata <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  
  
  # APPLY FILTER
  #----------------------------------------------------------------------------#
  filtered_taxa <- apply_filter_verb(metadata, long_abun, expressions, .by, "taxa") 

  filtered_taxa <- purrr::reduce(filtered_taxa, dplyr::semi_join, by = c("mgnet", "taxa_id")) 
  for(i in names(object)){
    filtered_i <- dplyr::filter(filtered_taxa, mgnet == i) %>% dplyr::pull("taxa_id")
    filtered_i <- filtered_i[order(which(taxa_id(object[[i]]) %in% filtered_i))]
    object[[i]] <- object[[i]][, filtered_i]
  }
  
  validObject(object)
  return(object)
})
