# CLASS
#------------------------------------------------------------------------------#
#' @title MetaGenomic NETwork List Class (mgnetList)
#' 
#' @description
#' An S4 class designed to hold a list of `mgnet` objects. This class ensures that all elements
#' within the list are instances of the `mgnet` class, providing a structured and type-safe way
#' to manage collections of metagenomic network data.
#'
#' @slot mgnets A list where each element is an instance of the `mgnet` class.
#'
#' @seealso \code{\link{mgnet}} for details on the `mgnet` class.
#' 
#' @name mgnetList-class
#' @rdname mgnetList-class
#' @exportClass mgnetList
setClass(
  "mgnetList",
  slots = c(mgnets = "list"),
  prototype = prototype(mgnets = list()),
  validity = function(object) {
    
    if (!all(sapply(object@mgnets, inherits, "mgnet"))) {
      return("All elements of `mgnets` must be instances of the `mgnet` class")
    }
    
    # Check if all elements in the list are named
    if (!all(nzchar(names(object@mgnets)))) {
      return("All elements in `mgnets` must be named.")
    }
    
    TRUE
    })

# CONSTRUCTOR
#------------------------------------------------------------------------------#
#' Create an mgnetList Object
#'
#' This function constructs an `mgnetList` object, which is designed to encapsulate and manage
#' a collection of `mgnet` objects in a structured list.
#'
#' @param mgnets An optional list of `mgnet` objects to initialize the `mgnetList`.
#'        Defaults to an empty list. Each element must be an instance of the `mgnet` class.
#' 
#' @return An `mgnets` object.
#' @export
#' @name mgnetList
mgnetList <- function(mgnets=list()) {
  if (!all(sapply(mgnets, inherits, "mgnet"))) {
    stop("All elements provided must be instances of the `mgnet` class")
  }
  new("mgnetList", mgnets=mgnets)
}