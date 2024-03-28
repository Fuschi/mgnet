# FILTER INFO SAMPLE
#------------------------------------------------------------------------------#
#' Filter Samples in mgnet Objects Based on info_sample Conditions
#'
#' @description
#' This function filters samples in an `mgnet` object or each `mgnet` object within
#' an `mgnetList` based on specified conditions applied to the `info_sample` data frame.
#' The conditions for filtering are based on expressions applied to columns within
#' the `info_sample` data frame.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Conditions to filter samples, which are passed unquoted and can
#' include dplyr-style filtering expressions. These conditions are applied to the
#' columns of the `info_sample` data frame, which is accessed as a tibble with
#' `sample_id` serving as the key column.
#'
#' @return An `mgnet` or `mgnetList` object with samples filtered according to the
#' specified conditions. The filtered object(s) retain only those samples that meet
#' the filtering criteria specified by the `...` arguments.
#'
#' @details The `filter_info_sample` function leverages the `info_sample` getter function
#' to access the metadata of samples as a tibble, where the row names are transformed
#' into a `sample_id` column. This enables seamless integration with dplyr's `filter`
#' function for specifying complex filtering conditions. The function then uses the
#' `sample_id` column to identify which samples should be retained in the filtered `mgnet`
#' or `mgnetList` object.
#'
#' @export
#' @importFrom dplyr filter
#' @importFrom rlang enquos
#' @name filter_info_sample
#' @aliases filter_info_sample,mgnet-method filter_info_sample,mgnetList-method
#' @seealso \code{\link[dplyr]{filter}} for details on filter conditions.
setGeneric("filter_info_sample", function(object, ...) {
  standardGeneric("filter_info_sample")
})

setMethod("filter_info_sample", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Filter the info_sample data frame based on the conditions
  filtered_info_sample <- dplyr::filter(info_sample(object, .fmt = "tbl"), !!!conditions)
  
  # Check if the filtered result is empty
  if(nrow(filtered_info_sample) == 0) {
    stop("No samples meet the specified filtering conditions. Please revise your criteria.")
  }
  
  # Extract the sample_id of the filtered samples
  filtered_sample_ids <- filtered_info_sample$sample_id
  
  # Subset the object using the mgnet extractor method based on filtered sample IDs
  subsetted_object <- object[filtered_sample_ids,]
  
  return(subsetted_object)
})


setMethod("filter_info_sample", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    
    # Filter the info_sample data frame based on the conditions
    filtered_info_sample <- dplyr::filter(info_sample(mgnet_obj, .fmt = "tbl"), !!!conditions)
    
    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_info_sample) == 0) {
      stop("No samples meet the specified filtering conditions in one or more mgnet objects. Please revise your criteria.")
    }
    
    # Extract the sample_id of the filtered samples
    filtered_sample_ids <- filtered_info_sample$sample_id
    
    # Subset the object using the mgnet extractor method based on filtered sample IDs
    subsetted_object <- mgnet_obj[filtered_sample_ids,]
    
    return(subsetted_object)
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(object)
})


# FILTER INFO TAXA
#------------------------------------------------------------------------------#
#' Filter Taxa in mgnet Objects Based on Taxonomic and Metadata Conditions
#'
#' @description
#' This function filters taxa in an `mgnet` object or each `mgnet` object within
#' an `mgnetList` based on specified conditions applied to the merged data from
#' the `info_taxa` and `lineage` slots. This enables filtering taxa using both
#' taxonomic classification and additional taxa metadata.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Conditions for filtering taxa, passed unquoted and allowing for
#' dplyr-style filtering expressions. These conditions apply to columns in the
#' merged dataset from `info_taxa` and `lineage`.
#'
#' @return An `mgnet` or `mgnetList` object with the taxa subset according to the
#' specified conditions. The structure of the object remains intact, but taxa that
#' do not meet the conditions are removed from the abundance data, lineage information,
#' and potentially network and community analysis results.
#'
#' @details The function first merges `info_taxa` data with `lineage` data to create
#' a comprehensive taxonomic dataset. It then filters this dataset based on the
#' provided conditions, identifying which taxa to retain. The resulting filtered
#' taxa list is used to subset the original `mgnet` or `mgnetList` object, thus
#' affecting related abundance and network data but not altering the sample metadata.
#'
#' @export
#' @importFrom dplyr filter inner_join
#' @importFrom methods validObject
#' @importFrom rlang enquos
#' @name filter_info_taxa
#' @aliases filter_info_taxa,mgnet-method filter_info_taxa,mgnetList-method
#' @seealso \code{\link[dplyr]{filter}} for details on filter expressions.
setGeneric("filter_info_taxa", function(object, ...) {standardGeneric("filter_info_taxa")})

setMethod("filter_info_taxa", "mgnet", function(object, ...) {
  
  if(length(object@info_taxa)==0 & length(object@lineage)==0){
    stop("slots info_taxa and lineage missing")
  }
  
  # Capture filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Merge 'info_taxa' and 'lineage' into a comprehensive dataset
  if(length(object@info_taxa)!=0 & length(object@lineage)!=0){
    
    merged_taxa <- inner_join(info_taxa(object, .fmt = "tbl"),
                              lineage(object, .fmt = "tbl"),
                              by = "taxa_id")
    
  } else if (length(object@info_taxa)==0 & length(object@lineage)!=0) {
    
    merged_taxa <- lineage(object, .fmt = "tbl")
    
  } else if (length(object@info_taxa)!=0 & length(object@lineage)==0) {
    
    merged_taxa <- info_taxa(object, .fmt = "tbl")
    
  }
  
  
  # Apply conditions to filter taxa
  filtered_taxa <- dplyr::filter(merged_taxa, !!!conditions)
  
  # Check if the filtered result is empty
  if(nrow(filtered_taxa) == 0) {
    stop("No taxa meet the specified filtering conditions. Please revise your criteria.")
  }
  
  # Extract taxa IDs that meet the filtering conditions
  filtered_taxa_ids <- filtered_taxa$taxa_id
  
  # Subset the mgnet object to retain only filtered taxa
  subsetted_object <- object[, filtered_taxa_ids]
  
  return(subsetted_object)
})

setMethod("filter_info_taxa", "mgnetList", function(object, ...) {
  # Capture filtering conditions as quosures for list application
  conditions <- rlang::enquos(...)
  
  # Apply filtering to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    filtered_taxa_obj <- filter_info_taxa(mgnet_obj, !!!conditions)
    
    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_taxa_obj) == 0) { # Adjust this condition as needed based on actual implementation
      stop("No taxa meet the specified filtering conditions in one or more mgnet objects. Please revise your criteria.")
    }
    
    return(filtered_taxa_obj)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  validObject(object)
  return(object)
})