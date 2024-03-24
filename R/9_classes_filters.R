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
    
    # Extract the sample_id of the filtered samples
    filtered_sample_ids <- filtered_info_sample$sample_id
    
    # Subset the object using the mgnet extractor method based on filtered sample IDs
    subsetted_object <- mgnet_obj[filtered_sample_ids,]
    
    return(subsetted_object)
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(object)
})