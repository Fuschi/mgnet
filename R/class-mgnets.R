#' @include class-mgnet.R
NULL

#' Class `mgnets`: a collection of `mgnet` objects
#'
#' @description
#' `mgnets` is an S4 container class designed to store and manage a named
#' collection of [`mgnet`] objects.
#'
#' It provides a consistent, type-safe interface for workflows involving
#' multiple metagenomic datasets, while preserving the integrity of each
#' contained `mgnet` object.
#'
#' Grouping variables stored at the `mgnets` level are interpreted at the
#' collection level. In this context, available grouping variables are computed
#' on the union of metadata columns across the contained `mgnet` objects,
#' assuming that downstream operations harmonize missing columns through
#' `dplyr::bind_rows()`, with absent columns filled with `NA`.
#'
#' @slot mgnets
#' A list of objects of class [`mgnet`]. If not empty, the list must be named
#' and names must be unique.
#'
#' @section Validity checks:
#' If the collection is not empty:
#' \itemize{
#'   \item all elements must inherit from class `mgnet`;
#'   \item all elements must have non-empty names;
#'   \item element names must be unique;
#'   \item each contained `mgnet` object must itself be valid;
#'   \item if the attribute `mgnet_groups` is present, it must be a unique,
#'   non-empty character vector whose values belong to the union of variables
#'   available across the collection.
#' }
#'
#' @seealso
#' [`mgnet()`] for the constructor of single `mgnet` objects.
#'
#' @name mgnets-class
#' @rdname mgnets-class
#' @exportClass mgnets
#' @importFrom methods setClass setValidity
setClass(
  "mgnets",
  slots = c(mgnets = "list"),
  prototype = prototype(mgnets = list()),
  validity = function(object) {
    if (length(object@mgnets) == 0L) return(TRUE)
    
    # All elements must be mgnet objects
    if (!all(vapply(object@mgnets, function(z) inherits(z, "mgnet"), logical(1)))) {
      cli::cli_abort(
        "All elements of {.cls mgnets} must be instances of the {.cls mgnet} class."
      )
    }
    
    # Names must exist and be non-empty
    nm <- names(object@mgnets)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort(
        "All elements of {.cls mgnets} must be named with non-empty strings."
      )
    }
    
    # Names must be unique
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort(
        "All element names in {.cls mgnets} must be unique."
      )
    }
    
    # Each inner mgnet must be valid
    lapply(object@mgnets, methods::validObject)
    
    # Consistency of the "mgnet_groups" attribute at collection level
    groups <- attr(object, "mgnet_groups", exact = TRUE)
    
    if (!is.null(groups)) {
      if (!is.character(groups)) {
        cli::cli_abort(
          c(
            "x" = 'Attribute {.field "mgnet_groups"} on {.cls mgnets} must be a character vector.',
            "i" = "Use {.fn group_mgnet}() to set grouping variables."
          )
        )
      }
      
      if (anyNA(groups) || any(!nzchar(groups))) {
        cli::cli_abort(
          'Grouping variables in {.field "mgnet_groups"} on {.cls mgnets} must be non-NA, non-empty strings.'
        )
      }
      
      if (anyDuplicated(groups) > 0L) {
        cli::cli_abort(
          'Grouping variables in {.field "mgnet_groups"} on {.cls mgnets} must be unique.'
        )
      }
      
      # Union of variables available across the collection;
      # downstream bind_rows() harmonizes missing columns with NA
      available <- meta_taxa_vars(object, .fmt = "unique")
      missing   <- setdiff(groups, available)
      
      if (length(missing)) {
        cli::cli_abort(
          c(
            "x" = 'Unknown grouping variable{?s} stored in attribute {.field "mgnet_groups"} on {.cls mgnets}: {.val {missing}}.',
            "v" = "Available variables across the collection are: {.val {available}}."
          )
        )
      }
    }
    
    TRUE
  }
)

#' Create an `mgnets` object
#'
#' @description
#' Construct an `mgnets` object from multiple [`mgnet`] objects, or from a
#' single named list of `mgnet` objects.
#'
#' This constructor performs basic input checks before creating the S4 object.
#' If the collection is not empty, all elements must be valid `mgnet` objects,
#' all names must be present and non-empty, and names must be unique.
#'
#' @param ... `mgnet` objects, or a single named list of `mgnet` objects.
#'
#' @return
#' An object of class `mgnets`.
#'
#' @seealso
#' [`mgnet()`] for creating individual `mgnet` objects.
#'
#' @importFrom methods new
#' @export
#' @name mgnets
mgnets <- function(...) {
  x <- list(...)
  
  # Single list argument
  if (length(x) == 1L && is.list(x[[1L]])) {
    x <- x[[1L]]
  }
  
  if (length(x) > 0L) {
    # Type check
    if (!all(vapply(x, function(z) inherits(z, "mgnet"), logical(1)))) {
      cli::cli_abort(
        "All inputs must be {.cls mgnet} objects or a list of {.cls mgnet}."
      )
    }
    
    # Names must be present and non-empty
    nm <- names(x)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort(
        "All mgnet objects must be named with non-empty strings."
      )
    }
    
    # Names must be unique
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort("Names must be unique.")
    }
  }
  
  methods::new("mgnets", mgnets = x)
}

#------------------------------------------------------------------------------#
# `mgnets` general methods
#------------------------------------------------------------------------------#

#' General methods for `mgnets` objects
#'
#' @description
#' Basic utilities for interacting with `mgnets` objects, including coercion to
#' a base list, retrieval of length, and access to element names.
#'
#' @section Methods:
#' \describe{
#'   \item{`as.list(x)`}{Convert an `mgnets` object to a plain list of `mgnet` objects.}
#'   \item{`as(x, "list")`}{Coerce an `mgnets` object to a base list.}
#'   \item{`length(x)`}{Return the number of contained `mgnet` objects.}
#'   \item{`names(x)`}{Return the names of the contained `mgnet` objects.}
#'   \item{`names(x) <- value`}{Replace the names of the contained `mgnet` objects.}
#' }
#'
#' @param x An object of class `mgnets`.
#' @param value A character vector used in the replacement form `names<-`.
#' @param ... Not used.
#'
#' @return
#' \itemize{
#'   \item `as.list()` and `as(., "list")` return a list of `mgnet` objects.
#'   \item `length()` returns an integer of length 1.
#'   \item `names()` returns a character vector.
#'   \item `names<-()` returns the updated `mgnets` object.
#' }
#'
#' @aliases
#' as.list.mgnets
#' coerce,mgnets,list-method
#' length,mgnets-method
#' names,mgnets-method
#' names<-,mgnets,character-method
#'
#' @importFrom methods setAs setMethod setReplaceMethod validObject
#' @name mgnets-methods
#' @rdname mgnets-methods
NULL

#' @rdname mgnets-methods
#' @export
as.list.mgnets <- function(x, ...) {
  x@mgnets
}

# No roxygen block needed
methods::setAs("mgnets", "list", function(from) from@mgnets)

#' @rdname mgnets-methods
#' @export
setMethod("length", "mgnets", function(x) {
  length(x@mgnets)
})

#' @rdname mgnets-methods
#' @export
setMethod("names", "mgnets", function(x) {
  names(x@mgnets)
})

#' @rdname mgnets-methods
#' @export
setReplaceMethod(
  "names",
  signature(x = "mgnets", value = "character"),
  function(x, value) {
    names(x@mgnets) <- value
    methods::validObject(x)
    x
  }
)

#------------------------------------------------------------------------------#
# List-like access for `mgnets`: $, $<-, [[, [[<-
#------------------------------------------------------------------------------#

#' Access or assign `mgnet` elements in an `mgnets` object
#'
#' @description
#' List-like accessors for `mgnets` objects using `$`, `$<-`, `[[`, and `[[<-`.
#'
#' These methods allow convenient extraction and replacement of individual
#' `mgnet` objects stored in the collection.
#'
#' @param x An object of class `mgnets`.
#' @param name A single character string giving the element name for `$` and `$<-`.
#' @param i An index or name for `[[` and `[[<-`.
#' @param j Ignored. Included only for signature compatibility.
#' @param value An object of class `mgnet` to assign.
#' @param ... Not used.
#'
#' @return
#' \itemize{
#'   \item `$` and `[[` return a single `mgnet` object.
#'   \item `$<-` and `[[<-` return the updated `mgnets` object.
#' }
#'
#' @aliases
#' $,mgnets-method
#' $<-,mgnets,mgnet-method
#' [[,mgnets-method
#' [[<-,mgnets,ANY,ANY,mgnet-method
#'
#' @name mgnets-access
#' @rdname mgnets-access
NULL

#' @rdname mgnets-access
#' @export
setMethod("$", "mgnets", function(x, name) {
  if (!name %in% names(x@mgnets)) {
    cli::cli_abort(
      "No {.cls mgnet} object named {.val {name}} found in the {.cls mgnets} object."
    )
  }
  x@mgnets[[name]]
})

#' @rdname mgnets-access
#' @export
methods::setReplaceMethod(
  "$",
  signature(x = "mgnets", value = "mgnet"),
  function(x, name, value) {
    if (!inherits(value, "mgnet")) {
      cli::cli_abort("Assigned value must be an {.cls mgnet} object.")
    }
    x@mgnets[[name]] <- value
    methods::validObject(x)
    x
  }
)

#' @rdname mgnets-access
#' @export
setMethod("[[", "mgnets", function(x, i, j, ...) {
  x@mgnets[[i]]
})

#' @rdname mgnets-access
#' @export
methods::setReplaceMethod(
  "[[",
  signature(x = "mgnets", i = "ANY", j = "ANY", value = "mgnet"),
  function(x, i, j, value) {
    if (!inherits(value, "mgnet")) {
      cli::cli_abort("Assigned value must be an {.cls mgnet} object.")
    }
    x@mgnets[[i]] <- value
    methods::validObject(x)
    x
  }
)