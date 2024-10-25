#------------------------------------------------------------------------------#
#' Initialize Sample Information
#'
#' This internal function initializes the sample information for an `mgnet` object.
#' If metadata is present in the object, it returns the metadata in tibble format.
#' If no metadata is available, it creates a tibble using the `sample_id` of the object.
#'
#' @param object An `mgnet` object.
#'
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @return A tibble containing either the sample metadata or just the `sample_id` column.
#' @keywords internal
initialize_meta <- function(object) {
  
  if(inherits(object, "mgnet")){
    
    if(miss_sample(object)) stop("Error: No sample available.")
    
    if (has_slot(object, "meta")) {
      return(meta(object, .fmt = "tbl"))
    } else {
      return(tibble::tibble(sample_id = sample_id(object)))
    }
    
  } else if(inherits(object, "mgnetList")){
    
    if(miss_sample(object, "any")) stop("Error: No sample available in at least one of the mgnet objects.")
    
    object %>%
      purrr::map(initialize_meta) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      return()
    
  } else {
    
    stop("Error: object could be only mgnet or mgnetList")
    
  }
}

#------------------------------------------------------------------------------#
#' Initialize Taxa Information
#'
#' This internal function initializes the taxa information for an `mgnet` object.
#' If metadata is present in the object, it returns the metadata in tibble format.
#' If no metadata is available, it creates a tibble using the `taxa_id` of the object.
#'
#' @param object An `mgnet` object.
#'
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @return A tibble containing either the sample metadata or just the `taxa_id` column.
#' @keywords internal
initialize_taxa <- function(object) {
  
  if(inherits(object, "mgnet")){
    
    if(miss_taxa(object)) stop("Error: No taxa available.")
    
    if (has_metataxa(object)) {
      return(taxa(object, .fmt = "tbl"))
    } else {
      return(tibble::tibble(taxa_id = taxa_id(object)))
    }
    
  } else if(inherits(object, "mgnetList")){
    
    if(miss_taxa(object, "any")) stop("Error: No taxa available in at least one of the mgnet objects.")
    
    object %>%
      purrr::map(initialize_taxa) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      return()
    
  } else {
    
    stop("Error: object could be only mgnet or mgnetList")
    
  }
}


#------------------------------------------------------------------------------#
#' Join and Convert Abundance Data to Long Format
#'
#' This internal function creates a grid of `sample_id` and `taxa_id` combinations 
#' and joins all required abundance variables in long format. It converts abundance data 
#' (e.g., `abun`, `rela`, `norm`) from wide to long format for easier manipulation.
#'
#' @param object An `mgnet` object.
#' @param needed_keys A list containing the abundance keys to be used (e.g., `abun`, `rela`, `norm`).
#'
#' @return A tibble with `sample_id`, `taxa_id`, and the abundance variables in long format.
#' @keywords internal
long_abundance_join <- function(object, needed_keys) {
  
  if(inherits(object, "mgnet")){
    
    if(miss_sample(object) || miss_taxa(object)) stop("Error: No sample or taxa available.")
    if(length(needed_keys) == 0) return(tibble::tibble())
  
      long_abun <- tidyr::expand_grid(
        sample_id = sample_id(object),
        taxa_id = taxa_id(object))
      
      for (abundance_key in needed_keys) {
        long_abun <- long_abun %>%
          dplyr::left_join(
            methods::slot(object, abundance_key) %>%
              tibble::as_tibble(rownames = "sample_id") %>%
              tidyr::pivot_longer(-sample_id,
                                  names_to = "taxa_id",
                                  values_to = abundance_key),
            by = c("sample_id", "taxa_id")
          )}
  
      return(long_abun)
    
  } else if (inherits(object, "mgnetList")){

    if(miss_sample(object, "any") || miss_taxa(object, "any")) stop("Error: No sample or taxa available in at least one of the mgnet objects.")
    
    long_abun <- purrr::map(object, long_abundance_join, needed_keys = needed_keys) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
    
    return(long_abun)
    
  } else {
    
    stop("Error: object could be only mgnet or mgnetList")
    
  }
}


#------------------------------------------------------------------------------#
#' Split and Arrange Sample Info by mgnet
#'
#' This internal function splits the mutated sample info by `mgnet`, arranges 
#' the samples to match their order in the corresponding `mgnet` object, and 
#' returns the result with sample IDs as row names.
#'
#' @param merged_meta A tibble containing the merged sample info with `mgnet` column.
#' @param object An `mgnetList` object to match the sample order.
#' 
#' @importFrom purrr imap
#' @importFrom dplyr arrange select
#' @importFrom tibble column_to_rownames
#'
#' @return A list of tibbles with sample IDs as row names, one for each `mgnet` object.
#' @keywords internal
split_arrange_merged_meta <- function(merged_meta, object) {
  
  splitted_meta <- merged_meta %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x, y) {
      dplyr::arrange(x, base::match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x) {
      x %>%
        dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  # Ensure the split list retains the original order of mgnet objects
  splitted_meta <- splitted_meta[names(object)]
  
  return(splitted_meta)
}


#------------------------------------------------------------------------------#
#' Split and Arrange Taxa Info by mgnet
#'
#' This internal function splits the mutated taxa info by `mgnet`, arranges 
#' the taxa to match their order in the corresponding `mgnet` object, and 
#' returns the result with taxa IDs as row names.
#'
#' @param merged_taxa A tibble containing the merged taxa info with `mgnet` column.
#' @param object An `mgnetList` object to match the sample order.
#' 
#' @importFrom purrr imap
#' @importFrom dplyr arrange select
#' @importFrom tibble column_to_rownames
#' @importFrom tidyselect any_of
#'
#' @return A list of tibbles with taxa IDs as row names, one for each `mgnet` object.
#' @keywords internal
split_arrange_merged_taxa <- function(merged_taxa, object) {
  
  splitted_taxa <- merged_taxa %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x, y) {
      dplyr::arrange(x, base::match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x) {
      x %>%
        dplyr::select(-tidyselect::any_of(c("mgnet", "comm_id"))) %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  # Ensure the split list retains the original order of mgnet objects
  splitted_taxa <- splitted_taxa[names(object)]
  
  return(splitted_taxa)
}
