#------------------------------------------------------------------------------#
#' Filter Sample Metadata in mgnet or mgnetList Objects
#'
#' This function applies specified filtering criteria exclusively to the metadata of samples within 
#' an `mgnetList` object. It allows for advanced querying based on metadata fields and supports
#' tidyverse-style syntax for manipulating these data.
#'
#' @param object An object of class \code{mgnetList} containing the data to be filtered.
#' @param ... Expressions that define the filtering criteria using variables available
#'            in the samples metadata of the object.
#' @param .by A character vector specifying the grouping variables for filtering, which defaults to
#'            \code{"sample_id"} for `mgnet` objects and \code{"mgnet", "sample_id"} for `mgnetList` objects.
#'            This parameter should not include 'taxa_id' as it focuses solely on sample metadata.
#'
#' @return Returns a modified \code{mgnetList} object containing only the samples that meet the specified 
#'         filtering criteria. The method operates on samples metadata within each `mgnet` object in the list.
#'
#' @details
#' The function uses several reserved keywords:
#'   - \code{mgnet}: Refers to each individual `mgnet` object within an `mgnetList`, used for handling 
#'     data across multiple datasets.
#'   - \code{sample_id}: Unique identifier for each taxon.
#'   - \code{comm_id}: Identifies community membership; used if community data is present.
#' Additional metadata columns within the taxa structure can also be utilized for detailed filtering.
#' This function facilitates complex filtering operations tailored to specific study requirements using
#' these reserved keywords.
#'
#' @export
#' @aliases filter_taxa_info,mgnetList-method
#' @importFrom dplyr filter group_by ungroup arrange left_join select distinct
#' @importFrom rlang enquos eval_tidy syms
#' @importFrom tibble as_tibble add_column
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map list_rbind imap
#' @importFrom methods slot
#' @importFrom magrittr %>%
setGeneric("filter_sample_info", function(object, ..., .by) {standardGeneric("filter_sample_info")})

setMethod("filter_sample_info", "mgnet", function(object, ..., .by) {
  
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
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  if(length(meta(object)) != 0){
    sample_info_filtered <- meta(object, .fmt = "tbl")
  } else {
    sample_info_filtered <- tibble::tibble(sample_id = sample_id(object))
  }
  
  # FILTER
  #----------------------------------------------------------------------------#
  sample_info_filtered <- sample_info_filtered %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(sample_id)
  
  return(object[sample_info_filtered, ])
  
})

#------------------------------------------------------------------------------#
setMethod("filter_sample_info", "mgnetList", function(object, ..., .by){
  
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
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_sample_filtered_merged <- object %>%
    purrr::map(\(x) {
      if(length(meta(x)) != 0){
        meta(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # FILTER
  #----------------------------------------------------------------------------#
  info_sample_filtered_merged <- info_sample_filtered_merged %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(mgnet, sample_id)
  
  # SPLIT FILTERD
  #----------------------------------------------------------------------------#
  object <- purrr::imap(object, \(x, name){
    filtered_samples_name <- info_sample_filtered_merged %>%
      dplyr::filter(mgnet == name) %>%
      dplyr::pull(sample_id) 
    return(object[[name]][filtered_samples_name, ])
  })
  
  return(object)
})


#------------------------------------------------------------------------------#
#' Filter Taxa Metadata in mgnet or mgnetList Objects
#'
#' This function applies specified filtering criteria exclusively to the metadata of taxa within 
#' an `mgnetList` object. It allows for advanced querying based on metadata fields and supports
#' tidyverse-style syntax for manipulating these data.
#'
#' @param object An object of class \code{mgnetList} containing the microbiome data to be filtered.
#' @param ... Expressions that define the filtering criteria using variables available
#'            in the taxa metadata of the object.
#' @param .by A character vector specifying the grouping variables for filtering, which defaults to
#'            \code{"sample_id"} for `mgnet` objects and \code{"mgnet", "sample_id"} for `mgnetList` objects.
#'            This parameter should not include 'taxa_id' as it focuses solely on sample metadata.
#'
#' @return Returns a modified \code{mgnetList} object containing only the taxa that meet the specified 
#'         filtering criteria. The method operates on taxa metadata within each `mgnet` object in the list.
#'
#' @details
#' The function uses several reserved keywords:
#'   - \code{mgnet}: Refers to each individual `mgnet` object within an `mgnetList`, used for handling 
#'     data across multiple datasets.
#'   - \code{taxa_id}: Unique identifier for each taxon.
#'   - \code{comm_id}: Identifies community membership; used if community data is present.
#' Additional metadata columns within the taxa structure can also be utilized for detailed filtering.
#' This function facilitates complex filtering operations tailored to specific study requirements using
#' these reserved keywords.
#'
#' @export
#' @aliases filter_taxa_info,mgnetList-method
#' @importFrom dplyr filter group_by ungroup arrange left_join select distinct
#' @importFrom rlang enquos eval_tidy syms
#' @importFrom tibble as_tibble add_column
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map list_rbind imap
#' @importFrom methods slot
#' @importFrom magrittr %>%
setGeneric("filter_taxa_info", function(object, ..., .by) {standardGeneric("filter_taxa_info")})

setMethod("filter_taxa_info", "mgnet", function(object, ..., .by) {
  
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
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  if(length(taxa(object)) != 0){
    taxa_info_filtered <- taxa(object, .fmt = "tbl")
  } else {
    taxa_info_filtered <- tibble::tibble(taxa_id = taxa_id(object))
  }
  
  # FILTER
  #----------------------------------------------------------------------------#
  taxa_info_filtered <- taxa_info_filtered %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(taxa_id)
  
  return(object[, taxa_info_filtered])
  
})


setMethod("filter_taxa_info", "mgnetList", function(object, ..., .by){
  
  # Ensure there are taxa to process
  if (any(ntaxa(object) == 0)) {
    stop("Error: No taxa available in at least one provided 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","taxa_id")
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_taxa_filtered_merged <- object %>%
    purrr::map(\(x) {
      if(length(taxa(x)) != 0){
        taxa(x, .fmt = "tbl")
      } else {
        tibble::tibble(taxa_id = taxa_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # FILTER
  #----------------------------------------------------------------------------#
  info_taxa_filtered_merged <- info_taxa_filtered_merged %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::filter(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup() %>%
    dplyr::select(mgnet, taxa_id)
  
  # SPLIT FILTERD
  #----------------------------------------------------------------------------#
  object <- purrr::imap(object, \(x, name){
    filtered_taxa_name <- info_taxa_filtered_merged %>%
      dplyr::filter(mgnet == name) %>%
      dplyr::pull(taxa_id) 
    return(object[[name]][, filtered_taxa_name])
  })
  
  return(object)
})