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


# Internal function to check the assigned elements for the setter methods of mgnetList
#------------------------------------------------------------------------------#

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
  
  if (!inherits(object, "mgnetList")) stop(sprintf("%s must be an mgnetList.", objectName))
  
  # Check if both inputs are lists and value is named
  if (!is.list(value)) errors <- c(errors,sprintf("%s must be a list.", valueName))
  if (length(object@mgnets) != length(value)) errors <- c(errors, sprintf("Lengths of %s and %s lists must be equal.", objectName, valueName))
  if (is.null(names(value))) errors <- c(errors, sprintf("%s list must have named elements.", valueName))
  if (any(duplicated(value))) errors <- c(erros, sprintf("%s cannot have duplicated names", valueName))
  
  # Check if names of both lists are present, even in different order
  if (!is.null(names(object@mgnets)) && (!setequal(names(object@mgnets), names(value)))) {
    errors <- c(errors, sprintf("The %s and %s lists must have the same names, but their order can differ.", objectName, valueName))
  }
  
  if(length(errors)!=0){
    errors <- paste("- ", errors, "\n", sep="")
    errors <- paste(errors, collapse="")
    errors <- paste("\nDEBUGGER of assign in",objectName,":\n", errors, sep="")
    stop(errors)
  }
  
  TRUE
}


#' Check Tibble in Assign Methods for mgnetList
#' 
#' @description Internal function to check if the tibble provided as a value in assign methods
#' is coherent with the `mgnetList` object. The value will replace all the taxa metadata in the `mgnetList`.
#' This function ensures that the provided value meets necessary conditions before assignment.
#' 
#' @param object An `mgnetList` object.
#' @param value A tibble or data frame that will be split into different tables, one for each taxa metadata.
#'        It is necessary that the columns `mgnet` and `sample_id`/`taxa_id` exist in the value.
#' @param sample_or_taxa Character that indicates if the assign tbl is for the sample or taxa metadata. Possible choices are "sample" or "taxa".
#' 
#' @return Logical; `TRUE` if all checks pass, otherwise an error is thrown.
#' @keywords internal
#' @export
is_assign_tbl <- function(object, value, sample_or_taxa) {
  # Capture object and value names
  objectName <- deparse(substitute(object))
  valueName <- deparse(substitute(value))
  sample_or_taxa <- match.arg(sample_or_taxa, c("sample", "taxa"))
  
  # Check if the object is an mgnetList
  if (!inherits(object, "mgnetList")) {
    stop(sprintf("%s must be an mgnetList.", objectName))
  }
  
  # Check if value is a tibble or data frame
  if (!inherits(value, c("tbl_df", "data.frame"))) {
    stop(sprintf("%s must be a tibble or data frame.", valueName))
  }
  
  # Check if required columns are present in value
  request_cols <- ifelse(sample_or_taxa=="sample", c("mgnet", "sample_id"), c("mgnet", "taxa_id"))
  value_cols <- colnames(value)
  missing_cols <- setdiff(request_cols, value_cols)
  
  if (length(missing_cols) > 0) {
    stop(sprintf("The following required columns are missing in %s: %s", 
                 valueName, paste(missing_cols, collapse = ", ")))
  }
  
  # Check if rownames are set
  if (has_rownames(value)) {
    stop(sprintf("In %s, the rownames must not be set.", valueName))
  }
  
  # Check if all unique elements in the 'mgnet' column are in the names of the object
  unique_mgnets <- unique(value$mgnet)
  object_mgnets <- names(object)
  missing_mgnets <- setdiff(unique_mgnets, object_mgnets)
  
  if (length(missing_mgnets) > 0) {
    stop(sprintf("The following 'mgnet' values in %s are not present in %s: %s", 
                 valueName, objectName, paste(missing_mgnets, collapse = ", ")))
  }
  
  # If all checks pass
  return(TRUE)
}




