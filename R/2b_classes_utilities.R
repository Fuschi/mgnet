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




