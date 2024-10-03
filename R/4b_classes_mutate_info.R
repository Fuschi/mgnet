#------------------------------------------------------------------------------#
#' Modify and Augment Sample Metadata in `mgnet` Objects
#'
#' This function dynamically manipulates the `meta` slot within `mgnet` objects by applying 
#' user-defined transformations to the sample metadata. It leverages the `tidyverse` tools 
#' such as `dplyr` to enable flexible data transformations.
#'
#' @param object An `mgnet` object. This function targets the `meta` slot, which contains 
#'        metadata for each sample in the object.
#' @param ... Dynamic expressions or functions to be applied to the metadata.
#'        These expressions can modify metadata fields in the `meta` slot.
#' @param .by A character vector specifying the columns to group the data by before applying 
#'        transformations. The default is `"sample_id"`. Grouping ensures that transformations 
#'        are contextually applied within each subgroup defined by `.by`. The use of `"taxa_id"` 
#'        as a grouping variable is strictly prohibited to maintain consistency.
#'
#' @return Returns the modified `mgnet` object with the transformed `meta` slot reflecting 
#'         the applied transformations. All other slots and structure within the `mgnet` 
#'         object remain unchanged.
#'
#' @export
#' @aliases mutate_sample_info,mgnet-method
#' @importFrom dplyr mutate group_by ungroup 
#' @importFrom rlang enquos syms eval_tidy quo_name
#' @importFrom tibble column_to_rownames tibble add_column
#' @importFrom purrr list_rbind map imap
setGeneric("mutate_sample_info", function(object, ..., .by) {standardGeneric("mutate_sample_info")})

setMethod("mutate_sample_info", "mgnet", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
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
    sample_info_mutated <- meta(object, .fmt = "tbl")
  } else {
    sample_info_mutated <- tibble::tibble(sample_id = sample_id(object))
  }
  
  # MUTATE
  #----------------------------------------------------------------------------#
  sample_info_mutated <- sample_info_mutated %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup()
  
  meta(object) <- tibble::column_to_rownames(sample_info_mutated, "sample_id")
  return(object)
  
})

#------------------------------------------------------------------------------#
setMethod("mutate_sample_info", "mgnetList", function(object, ..., .by){
  
  # Ensure there are samples to process
  if (any(nsample(object) == 0)) {
    stop("Error: No samples available in at least one provided 'mgnet' object.")
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
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","sample_id")
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_sample_mutated_merged <- object %>%
    purrr::map(\(x) {
      if(length(meta(x)) != 0){
        meta(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # MUTATE
  #----------------------------------------------------------------------------#
  info_sample_mutated_merged <- info_sample_mutated_merged %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup()
  
  # SPLIT MUTATED
  #----------------------------------------------------------------------------#
  info_sample_mutated_splitted <- info_sample_mutated_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  meta(object) <- info_sample_mutated_splitted
  return(object)
})


#------------------------------------------------------------------------------#
#' Modify and Augment Taxa Metadata in `mgnet` Objects
#'
#' This function dynamically manipulates the `taxa` slot within `mgnet` objects by applying 
#' user-defined transformations to the taxa metadata. It leverages the `tidyverse` tools 
#' such as `dplyr` to enable flexible data transformations.
#'
#' @param object An `mgnet` object. This function targets the `taxa` slot, which contains 
#'        metadata for each taxon in the object.
#' @param ... Dynamic expressions or functions to be applied to the metadata.
#'        These expressions can modify metadata fields in the `taxa` slot.
#' @param .by A character vector specifying the columns to group the data by before applying 
#'        transformations. The default is `"taxa_id"`. Grouping ensures that transformations 
#'        are contextually applied within each subgroup defined by `.by`. The use of `"sample_id"` 
#'        as a grouping variable is strictly prohibited to maintain consistency.
#'
#' @return Returns the modified `mgnet` object with the transformed `meta` slot reflecting 
#'         the applied transformations. All other slots and structure within the `mgnet` 
#'         object remain unchanged.
#'
#' @export
#' @aliases mutate_sample_info,mgnet-method
#' @importFrom dplyr mutate group_by ungroup 
#' @importFrom rlang enquos syms eval_tidy quo_name
#' @importFrom tibble column_to_rownames tibble add_column
#' @importFrom purrr list_rbind map imap
setGeneric("mutate_taxa_info", function(object, ..., .by) {standardGeneric("mutate_taxa_info")})

setMethod("mutate_taxa_info", "mgnet", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (ntaxa(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
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
    taxa_info_mutated <- taxa(object, .fmt = "tbl")
  } else {
    taxa_info_mutated <- tibble::tibble(taxa_id = taxa_id(object))
  }
  
  # MUTATE
  #----------------------------------------------------------------------------#
  taxa_info_mutated <- taxa_info_mutated %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup()
  
  taxa(object) <- tibble::column_to_rownames(taxa_info_mutated, "taxa_id")
  return(object)
  
})

#------------------------------------------------------------------------------#
setMethod("mutate_taxa_info", "mgnetList", function(object, ..., .by){
  
  # Ensure there are taxa to process
  if (any(ntaxa(object) == 0)) {
    stop("Error: No taxa available in at least one provided 'mgnet' object.")
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
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","taxa_id")
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_taxa_mutated_merged <- object %>%
    purrr::map(\(x) {
      if(length(taxa(x)) != 0){
        taxa(x, .fmt = "tbl")
      } else {
        tibble::tibble(taxa_id = taxa_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # MUTATE
  #----------------------------------------------------------------------------#
  info_taxa_mutated_merged <- info_taxa_mutated_merged %>%
    dplyr::group_by(!!!rlang::syms(.by)) %>%
    dplyr::mutate(!!!rlang::eval_tidy(expressions)) %>%
    dplyr::ungroup()
  
  # SPLIT MUTATED
  #----------------------------------------------------------------------------#
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