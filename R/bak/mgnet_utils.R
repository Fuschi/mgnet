#' Check if an object is an `mgnet` Instance
#'
#' This function checks whether a given object is an instance of the `mgnet` class.
#'
#' @param object The object to check.
#' @return `TRUE` if the object is an instance of the `mgnet` class, `FALSE` otherwise.
#' @examples
#' my_object <- mgnet(abundance = matrix(1:4, ncol = 2))
#' is_mgnet(my_object) # Returns TRUE
#' is_mgnet(list()) # Returns FALSE
#' @export
is_mgnet <- function(object) {inherits(object, "mgnet")}

#' Ensure Object is an `mgnet` Instance
#'
#' Checks if the provided object is an instance of the `mgnet` class
#' and stops execution with an error message if not.
#'
#' @param object The object to check.
#' @return Invisible `TRUE` if the object is an instance of the `mgnet` class.
#' @examples
#' my_object <- new("mgnet", abundance = matrix(1:4, ncol = 2))
#' ensure_mgnet(my_object) # Does nothing (object is of class 'mgnet')
#' ensure_mgnet(list()) # Throws an error
#' @export
ensure_mgnet <- function(object) {
  if (!is_mgnet(object)) {
    stop(sprintf("'%s' is not of class 'mgnet'.", deparse(substitute(object))))
  }
  invisible(TRUE)
}
