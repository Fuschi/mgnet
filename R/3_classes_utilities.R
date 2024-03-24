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


#' Check list in assign methods for mgnetList
#' 
#' @description Internal function to check if the list provided as a value in assign methods
#' is coherent with the mgnetList object.
#' @param object mgnetList object.
#' @param value List intended to be assigned to the mgnetList.
#' @keywords internal
#' @export
are_list_assign <- function(object, value) {
  
  objectName <- deparse(substitute(object))
  valueName <- deparse(substitute(value))
  errors <- character()
  
  # Check if both inputs are lists and value is named
  if (!is.list(value)) errors <- c(errors,sprintf("%s must be a list.", valueName))
  if (!inherits(object, "mgnetList")) errors <- c(errors, sprintf("%s must be an mgnetList.", objectName))
  if (length(object@mgnets) != length(value)) errors <- c(errors, sprintf("Lengths of %s and %s lists must be equal.", objectName, valueName))
  if (is.null(names(value))) errors <- c(errors, sprintf("%s list must have named elements.", valueName))
  
  # Check names alignment if both have names
  if (!is.null(names(object@mgnets)) && any(names(object@mgnets) != names(value))) {
    errors <- c(errors, sprintf("Elements of %s and %s must have the same names in the same order.", objectName, valueName))
  }
  
  if(length(errors)!=0){
    errors <- paste("- ", errors, "\n", sep="")
    errors <- paste(errors, collapse="")
    errors <- paste("\nDEBUGGER of assign in",objectName,":\n", errors, sep="")
    stop(errors)
  }
  
  TRUE
}


# UPDATE SAMPLE SUM (Information Necessary To Calculates The Relative Abundances)
#------------------------------------------------------------------------------#
#' Update Sample Sum in `mgnet` Objects
#'
#' This function updates the `sample_sum` field in the `info_sample` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`. The `sample_sum` represents the total count/abundance
#' for each sample, crucial for calculating relative abundances. For `mgnetList` objects, the function
#' iteratively applies this process to each contained `mgnet` object.
#'
#' @param object An `mgnet` or `mgnetList` object, or a list of `mgnet` objects. When the object is
#'        of class `mgnet`, `sampleSums` can optionally be provided to manually set the sample sums.
#'        For `mgnetList` objects and when `sampleSums` is missing, sample sums are automatically
#'        calculated based on the abundance data.
#' @param sampleSums A numeric vector of sample sums to be explicitly set for an `mgnet` object.
#'        This parameter is ignored for `mgnetList` objects and when calculating sample sums
#'        automatically. Optional for `mgnet` objects.
#'
#' @details For `mgnet` objects, if `sampleSums` is provided, these values are used to update the
#'          `sample_sum` in the `info_sample` slot. Otherwise, sample sums are automatically calculated
#'          from the abundance data. For `mgnetList` objects, sample sums are always automatically
#'          calculated for each contained `mgnet` object based on their abundance data.
#'
#' @return The modified `mgnet` or `mgnetList` object with updated `sample_sum` in the `info_sample` slot.
#'
#' @export
#' @name update_sample_sum
#' @aliases update_sample_sum,mgnet-method update_sample_sum,mgnetList-method
setGeneric("update_sample_sum", function(object, sampleSums) standardGeneric("update_sample_sum"))

#' @rdname update_sample_sum
setMethod("update_sample_sum", c("mgnet","missing"), function(object) {
  if(length(object@abundance) == 0) stop("abundance slot cannot be empty.")
  
  # Calculate sample sums and update or create the 'sample_sum' column
  sampleSums <- rowSums(object@abundance, na.rm = TRUE)
  
  if(length(object@info_sample)==0){
    info_sample <- data.frame("sample_sum"=sampleSums)
    rownames(info_sample) <- rownames(object@abundance)
  } else if( length(info_sample)!=0 && "sample_sum"%in%colnames(info_sample) ){
    message("The `sample_sum` column in `info_sample` slot has been updated.")
    info_sample$sample_sum <- sampleSums
  } else {
    info_sample$sample_sum <- sampleSums
  }
  
  validObject(object)
  return(object)
})

#' @rdname update_sample_sum
setMethod("update_sample_sum", c("mgnet","ANY"), function(object, sampleSums) {
  
  if(length(object@abundance)==0) stop("abundance slot cannot be empty.")
  if(!is.vector(sampleSums) || !is.numeric(sampleSums)) stop("sampleSums must be a numeric vector.")
  if(any(sampleSums<=0)) stop("sampleSums elements must be all > 0.")
  if(length(sampleSums)!=nrow(object@abundance)) stop("sampleSums length not match the number of samples (rows) of abundance.")
  
  if(length(object@info_sample)==0){
    info_sample <- data.frame("sample_sum"=sampleSums)
    rownames(info_sample) <- rownames(object@abundance)
  } else if( length(info_sample)!=0 && "sample_sum"%in%colnames(info_sample) ){
    message("The `sample_sum` column in `info_sample` slot has been updated.")
    info_sample$sample_sum <- sampleSums
  } else {
    info_sample$sample_sum <- sampleSums
  }
  
  validObject(object)
  return(object)
})

#' @rdname update_sample_sum
setMethod("update_sample_sum", c("mgnetList","missing"), function(object) {
  object@mgnets <- sapply(object@mgnets, update_sample_sum, simplify = F, USE.NAMES = T)
  validObject(object)
  return(object)
})

#' @rdname update_sample_sum
setMethod("update_sample_sum", c("mgnetList","list"), function(object, sampleSums) {
  are_list_assign(object, sampleSums)
  object@mgnets <- sapply(object@mgnets, update_sample_sum, simplify = F, USE.NAMES = T)
  validObject(object)
  return(object)
})
# END UPDATE SAMPLE SUM
#------------------------------------------------------------------------------#