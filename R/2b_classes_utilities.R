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


#' Convert an mgnetList Object to a List
#'
#' This method allows for the conversion of an `mgnetList` object into a list, 
#' facilitating the application of list-based operations and functions. Each element 
#' of the list represents an `mgnet` object contained within the `mgnetList`.
#'
#' @param x An `mgnetList` object.
#' @param ... Additional arguments affecting the conversion (currently unused).
#'
#' @return A list where each element is an `mgnet` object previously contained 
#'         within the `mgnetList`. The names of the list elements correspond to 
#'         the names of the `mgnet` objects in the `mgnetList`, if they are named.
#'
#' @importFrom methods slot
#' @export
#' @seealso \link{mgnetList}, for details on the `mgnetList` class.
#' @name as.list-mgnetList
#' @aliases as.list,mgnetList-method
setMethod("as.list", "mgnetList", function(x, ...){
  
  slot(x, "mgnets")
})


#' Length of an mgnetList
#'
#' @description
#' Returns the number of `mgnet` objects contained within an `mgnetList` object.
#'
#' @param x An `mgnetList` object.
#' @return Integer value representing the number of `mgnet` objects in the `mgnetList`.
#' @export
setMethod("length", "mgnetList", function(x) {
  length(x@mgnets)
})


#' Names of mgnet Objects in an mgnetList
#'
#' @description
#' Retrieves the names of `mgnet` objects contained within an `mgnetList`.
#' Names provide a convenient way to reference and manage individual `mgnet` objects.
#'
#' @param x An `mgnetList` object.
#' @return A character vector of names of the `mgnet` objects.
#'         
#' @export
setMethod("names", "mgnetList", function(x) {
  names(x@mgnets)
})

#' Set Names of mgnet Objects in an mgnetList
#'
#' @description
#' Sets the names of `mgnet` objects contained within an `mgnetList`.
#' Names provide a convenient way to reference and manage individual `mgnet` objects.
#'
#' @param x An `mgnetList` object.
#' @param value A character vector representing the new names to be assigned to the `mgnet` objects.
#' @return The modified `mgnetList` object with updated names.
#'         
#' @importFrom methods validObject
#' @export
setMethod("names<-", "mgnetList", function(x, value) {
  names(x@mgnets) <- value
  validObject(x)
  x
})




