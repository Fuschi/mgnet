#' MetaGenomic NETwork List (mgnetList) Class
#'
#' @description
#' mgnetList is an S4 class designed to hold and manage a collection of mgnet objects,
#' providing a structured way to handle multiple metagenomic networks. It ensures that all
#' elements within the list adhere to the mgnet class structure, making it useful for
#' analyses involving multiple metagenomic datasets. This class allows for consistent and
#' type-safe manipulation of multiple mgnet objects.
#'
#' @slot mgnets A list where each element is an instance of the mgnet class.
#' If the list is not empty, each mgnet object within the list must be named to facilitate
#' easy access and management.
#'
#' @section Validity Checks:
#' The class includes validity checks to ensure that, if not empty, all elements within the
#' mgnets slot are instances of the mgnet class and are named. These checks help maintain
#' integrity and usability, preventing errors during downstream analyses.
#'
#' @seealso \code{\link[=mgnet]{mgnet}} for details on the mgnet class and its functionalities.
#'
#' @name mgnetList-class
#' @rdname mgnetList-class
#' @exportClass mgnetList
setClass(
  "mgnetList",
  slots = c(mgnets = "list"),
  prototype = prototype(mgnets = list()),
  validity = function(object) {
    
    if(length(object@mgnets)!=0){
      
      if (!all(sapply(object@mgnets, inherits, "mgnet"))) {
        return("All elements of `mgnets` must be instances of the `mgnet` class")
      }
      
      
      # Check if all elements in the list are named
      if(is.null(names(object@mgnets))){
        return("All elements in `mgnets` must be named.")
      } else if (any(nzchar(names(object@mgnets)))) {
        return("All elements in `mgnets` must be named.")
      }
      
    }
      
    TRUE
  })


# CONSTRUCTOR
#------------------------------------------------------------------------------#
#' MetaGenomic NETwork List (mgnetList) Class
#'
#' @description
#' mgnetList is an S4 class designed to hold and manage a collection of mgnet objects,
#' providing a structured way to handle multiple metagenomic networks. It ensures that all
#' elements within the list adhere to the mgnet class structure, making it useful for
#' analyses involving multiple metagenomic datasets. This class allows for consistent and
#' type-safe manipulation of multiple mgnet objects.
#'
#' @slot mgnets A list where each element is an instance of the mgnet class.
#' If the list is not empty, each mgnet object within the list must be named to facilitate
#' easy access and management.
#'
#' @section Validity Checks:
#' The class includes validity checks to ensure that, if not empty, all elements within the
#' mgnets slot are instances of the mgnet class and are named. These checks help maintain
#' integrity and usability, preventing errors during downstream analyses.
#'
#' @seealso \code{\link[=mgnet]{mgnet}} for details on the mgnet class and its functionalities.
#'
#' @name mgnetList-class
#' @rdname mgnetList-class
#' @exportClass mgnetList
setClass(
  "mgnetList",
  slots = c(mgnets = "list"),
  prototype = prototype(mgnets = list()),
  validity = function(object) {
    
    # Check if list is not empty
    if(length(object@mgnets)!=0){
      
      # Ensure all elements are instances of mgnet
      if (!all(sapply(object@mgnets, inherits, "mgnet"))) {
        return("\nAll elements of mgnets must be instances of the mgnet class.")
      }
      # Ensure all elements are named
      if(is.null(names(object@mgnets)) || any(!nzchar(names(object@mgnets)))) {
        return("\nAll elements in mgnets must be named.")
      }
      # Ensure all names are unique
      if(length(unique(names(object@mgnets))) != length(object@mgnets)) {
        return("\nAll elements names must be unique.")
      }
      
    }
    TRUE
  }
)

#CONSTRUCTOR
#------------------------------------------------------------------------------#
#' Create an mgnetList Object
#'
#' @description
#' Constructs an mgnetList object to encapsulate and manage a collection of mgnet objects.
#' This function ensures that all provided mgnet objects are validated and organized into a
#' structured list, facilitating easy access to multiple metagenomic networks. The constructor
#' supports the initialization of an empty list or a list populated with mgnet objects.
#'
#' @param ... Optional mgnet objects to be included in the mgnetList.
#' Providing named mgnet objects enhances manageability and clarity. The list can be
#' initialized as empty.
#'
#' @return An mgnetList object containing any provided mgnet objects. If not empty,
#' each mgnet object must be named, ensuring structured and accessible management
#' of multiple metagenomic networks.
#'
#' @export
#' @name mgnetList
mgnetList <- function(...) {
  mgnets <- list(...)
  new("mgnetList", mgnets=mgnets)
}


#' Convert an mgnetList Object to a List
#'
#' This function converts an \code{mgnetList} object to a list by extracting the individual
#' \code{mgnet} objects contained within the list.
#'
#' @param x An \code{mgnetList} object to be converted to a list.
#' @param ... Additional arguments (not currently used).
#'
#' @return A list containing the individual \code{mgnet} objects from the \code{mgnetList}.
#'
#' @seealso \code{\link{as.list}} for the generic method to convert objects to lists.
#' @export
setMethod("as.list", signature = "mgnetList", function(x, ...) {
  x@mgnets
})
