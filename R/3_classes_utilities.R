#' Check if an object is an `mgnet` Instance
#'
#' This function checks whether a given object is an instance of the `mgnet` class.
#'
#' @param object The object to check.
#' @return `TRUE` if the object is an instance of the `mgnet` class, `FALSE` otherwise.
#' @export
is_mgnet <- function(object) {inherits(object, "mgnet")}

#' Ensure Object is an `mgnet` Instance
#'
#' Checks if the provided object is an instance of the `mgnet` class
#' and stops execution with an error message if not.
#'
#' @param object The object to check.
#' @return Invisible `TRUE` if the object is an instance of the `mgnet` class.
#' @export
ensure_mgnet <- function(object) {
  if (!is_mgnet(object)) {
    stop(sprintf("'%s' is not of class 'mgnet'.", deparse(substitute(object))))
  }
  invisible(TRUE)
}

#' Check if an Object is an `mgnetList` Instance
#'
#' This function checks whether a given object is an instance of the `mgnetList` class.
#'
#' @param object The object to check.
#' @return `TRUE` if the object is an instance of the `mgnetList` class, `FALSE` otherwise.
#' @export
is_mgnetList <- function(object) {
  inherits(object, "mgnetList")
}

#' Ensure Object is an `mgnetList` Instance
#'
#' Checks if the provided object is an instance of the `mgnetList` class
#' and stops execution with an error message if not.
#'
#' @param object The object to check.
#' @return Invisible `TRUE` if the object is an instance of the `mgnetList` class.
#' @export
ensure_mgnetList <- function(object) {
  if (!is_mgnetList(object)) {
    stop(sprintf("'%s' is not of class 'mgnetList'.", deparse(substitute(object))))
  }
  invisible(TRUE)
}

# UPDATE SAMPLE SUM (Information Necessary To Calculates The Relative Abundances)
#------------------------------------------------------------------------------#
#' Update Sample Sum in `mgnet` Objects
#'
#' This function updates the `sample_sum` field in the `info_sample` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`. The `sample_sum` represents the total count/abundance
#' for each sample, which is crucial for calculating relative abundances. Note that `sample_sum`,
#' along with `taxa_id` and `sample_id`, is a reserved keyword within the `mgnet` class and should
#' not be used as a column name in the `info_sample` or any other slot externally.
#'
#' @param object An object of class `mgnet` or a list of `mgnet` objects contained within an `mgnetList`.
#'
#' @details For `mgnet` objects, the function calculates the row sums of the `abundance` slot, storing
#' these sums in the `info_sample` slot under the column `sample_sum`. If `info_sample` does not exist,
#' it initializes this slot as a data.frame with the `sample_sum` column. For `mgnetList` objects,
#' it iteratively applies this process to each contained `mgnet` object.
#'
#' @return The modified `mgnet` or `mgnetList` object with updated `sample_sum` in the `info_sample` slot.
#'
#' @export
#' @name update_sample_sum
#' @aliases update_sample_sum,mgnet-method update_sample_sum,mgnetList-method
setGeneric("update_sample_sum", function(object) standardGeneric("update_sample_sum"))

setMethod("update_sample_sum", "mgnet", function(object) {
  if(length(object@abundance) == 0) stop("abundance slot cannot be empty.")
  
  # Calculate sample sums and update or create the 'sample_sum' column
  sampleSums <- rowSums(object@abundance, na.rm = TRUE)
  
  if(length(object@info_sample)==0){
    info_sample <- data.frame("sample_sum"=sampleSums)
  } else if( length(info_sample)!=0 & !("sample_sum"%in%colnames(info_sample))){
    info_sample$sample_sum <- sampleSums
  }
  
  message("The 'sample_sum' column in 'info_sample' slot has been updated.")
  return(object)
})

setMethod("update_sample_sum", "mgnetList", function(object) {
  # Ensure the object contains mgnet objects
  if(!all(sapply(object@mgnetObjects, is_mgnet))) {
    stop("All elements in the list must be 'mgnet' objects.")
  }
  
  # Apply update_sample_sum to each mgnet object in the list
  object@mgnetObjects <- lapply(object@mgnetObjects, update_sample_sum)
  
  return(object)
})
# END UPDATE SAMPLE SUM
#------------------------------------------------------------------------------#