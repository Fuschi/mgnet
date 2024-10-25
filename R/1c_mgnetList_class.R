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
    errors <- character()
    mgnet_errors <- character()
    
    # Global mgnetList checks
    if (length(object@mgnets) != 0) {
      # Ensure all elements are instances of mgnet
      if (!all(sapply(object@mgnets, inherits, "mgnet"))) {
        errors <- c(errors, "All elements of `mgnets` must be instances of the `mgnet` class.")
      }
      
      # Ensure all elements are named
      if (is.null(names(object@mgnets)) || any(!nzchar(names(object@mgnets)))) {
        errors <- c(errors, "All elements in `mgnets` must be named.")
      }
      
      # Ensure names are unique
      if (length(unique(names(object@mgnets))) != length(object@mgnets)) {
        errors <- c(errors, "All element names must be unique.")
      }
      
      # Validate each mgnet object individually
      for (i in seq_along(object@mgnets)) {
        mgnet_obj <- object@mgnets[[i]]
        mgnet_name <- names(object@mgnets)[i]
        
        # Check if the individual mgnet object passes validation
        mgnet_validity <- tryCatch(
          {
            validObject(mgnet_obj, test = TRUE)
            TRUE  # Return TRUE if validObject passes without error
          },
          error = function(e) {
            e$message  # Capture error message if validity fails
          }
        )
        
        # If validation fails for mgnet, accumulate the error
        if (mgnet_validity != TRUE) {
          if (is.null(mgnet_name) || !nzchar(mgnet_name)) mgnet_name <- "NA"
          mgnet_errors <- c(mgnet_errors, sprintf("==== '%s' error report ====", mgnet_name), mgnet_validity)
        }
      }
      
      # Combine global errors and mgnet-specific errors
      all_errors <- character()
      if (length(errors) > 0) {
        errors <- lapply(errors, function(x) paste("- ", x, sep = ""))
        all_errors <- c("\n==== mgnetList error report ====", paste(errors, collapse = "\n"))
      }
      if (length(mgnet_errors) > 0) {
        all_errors <- c(all_errors, mgnet_errors)
      }
      
      # Return the combined error messages if any errors were found
      if (length(all_errors) > 0) {
        return(paste(all_errors, collapse = "\n"))
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



#------------------------------------------------------------------------------#
#' Access or Assign Elements in an mgnetList Object
#'
#' These methods allow accessing or assigning individual `mgnet` objects within an `mgnetList` 
#' using the `$` and `$<-` operators, providing list-like behavior for the `mgnetList` class.
#'
#' @section Usage:
#' \preformatted{
#' mgl$name
#' mgl$name <- value
#' }
#'
#' @param x An `mgnetList` object.
#' @param name A character string representing the name of the `mgnet` object to extract or assign.
#' @param value A `mgnet` object to be assigned to the `mgnetList`.
#' 
#' @details
#' \code{mgl$name} extracts the `mgnet` object from the `mgnetList` named \code{name}. If no 
#' such object exists, an error is thrown.
#' 
#' \code{mgl$name <- value} assigns a new `mgnet` object to the element named \code{name} in 
#' the `mgnetList`. If the value is not a valid `mgnet` object, an error is thrown. After the 
#' assignment, the validity of the entire `mgnetList` object is checked to ensure it adheres 
#' to the defined constraints.
#'
#' @return
#' - For `$`, the method returns the `mgnet` object corresponding to \code{name}.
#' - For `$<-`, the method assigns a new `mgnet` object to the list and returns the modified `mgnetList`.
#'
#' @examples
#' data(mgl, package = "mgnet")
#' 
#' # Access a specific mgnet object
#' mgnet_obj <- mgl$`69-001`
#' 
#' # Assign a new mgnet object to the mgnetList
#' mgl$new_sample <- mgnet()
#' 
#' @name mgnetList-access
#' @rdname mgnetList-access
#' @aliases $,mgnetList-method $<-,mgnetList-method
#' @export
setMethod("$", "mgnetList", function(x, name) {
  if (!name %in% names(x@mgnets)) {
    stop(sprintf("No mgnet object named '%s' found.", name))
  }
  x@mgnets[[name]]
})

#' @rdname mgnetList-access
#' @export
setMethod("$<-", "mgnetList", function(x, name, value) {
  if (!inherits(value, "mgnet")) {
    stop("Assigned value must be an `mgnet` object.")
  }
  
  x@mgnets[[name]] <- value
  
  # Re-validate the object after assignment
  validObject(x)
  
  x
})
