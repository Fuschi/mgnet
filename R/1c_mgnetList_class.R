# Define mgnetList class
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
    
    if (length(object@mgnets) != 0) {
      # Ensure all elements are instances of mgnet
      if (!all(sapply(object@mgnets, inherits, "mgnet"))) {
        return("All elements of `mgnets` must be instances of the `mgnet` class.")
      }
      
      # Ensure all elements are named
      if (is.null(names(object@mgnets)) || any(!nzchar(names(object@mgnets)))) {
        return("All elements in `mgnets` must be named.")
      }
      
      # Ensure all names are unique
      if (length(unique(names(object@mgnets))) != length(object@mgnets)) {
        return("All element names must be unique.")
      }
    }
    
    TRUE
  }
)

# Define constructor for mgnetList
#------------------------------------------------------------------------------#
#' Create an mgnetList Object
#'
#' @description
#' Constructs an mgnetList object to encapsulate and manage a collection of mgnet objects.
#' This function ensures that all provided mgnet objects are validated and organized into a
#' structured list, facilitating easy access to multiple metagenomic networks. The constructor
#' supports the initialization of an empty list, individual mgnet objects, or a list of mgnet objects.
#'
#' @param ... Optional mgnet objects or a list of mgnet objects to be included in the mgnetList.
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
  
  # Check if the first argument is a list of mgnet objects
  if (length(mgnets) == 1 && is.list(mgnets[[1]]) && all(sapply(mgnets[[1]], inherits, "mgnet"))) {
    mgnets <- mgnets[[1]]
  }
  
  # Ensure all provided mgnet objects are named
  if (length(mgnets) != 0 && (is.null(names(mgnets)) || any(!nzchar(names(mgnets))))) {
    stop("All mgnet objects must be named.")
  }
  
  new("mgnetList", mgnets = mgnets)
}

#' Convert an mgnetList Object to a List
#'
#' This function converts an `mgnetList` object to a list by extracting the individual
#' `mgnet` objects contained within the list.
#'
#' @param x An `mgnetList` object to be converted to a list.
#' @param ... Additional arguments (not currently used).
#' @return A list containing the individual `mgnet` objects from the `mgnetList`.
#' @export
as.list.mgnetList <- function(x, ...) {
  x@mgnets
}


# Define coercion from mgnetList to list
#------------------------------------------------------------------------------#
#' @export
setAs("mgnetList", "list", function(from) from@mgnets)
