#' @include class-mgnet.R
NULL

#' Class `mgnets` – a collection of `mgnet` objects
#'
#' @description
#' `mgnets` is an S4 class to hold and manage a collection of `mgnet` objects,
#' enabling multi-dataset analyses with a consistent, type-safe interface.
#'
#' @slot mgnets A list whose elements are instances of class `mgnet`.
#' If not empty, the list must be named and names must be unique.
#'
#' @section Validity checks:
#' If not empty:
#' - all elements must inherit from class `mgnet`;
#' - all elements must have non-empty names;
#' - names must be unique.
#'
#' @seealso [mgnet()] for details about the single-mgnet class used as elements.
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
    
    # all elements are mgnet
    if (!all(vapply(object@mgnets, function(z) inherits(z, "mgnet"), logical(1)))) {
      cli::cli_abort("All elements of {.cls mgnets} must be instances of the {.cls mgnet} class.")
    }
    # names exist and are non-empty
    nm <- names(object@mgnets)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort("All elements of {.cls mgnets} must be named with non-empty strings.")
    }
    # names are unique
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort("All element names in {.cls mgnets} must be unique.")
    }
    
    # 1) Ogni mgnet interno deve essere valido (inclusi group/link checks)
    lapply(object@mgnets, methods::validObject)
    
    # 2) Coerenza dell'attributo "mgnet_groups" a livello mgnets
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
          'Grouping variables in "mgnet_groups" on {.cls mgnets} must be non-NA, non-empty strings.'
        )
      }
      if (anyDuplicated(groups) > 0L) {
        cli::cli_abort(
          'Grouping variables in "mgnet_groups" on {.cls mgnets} must be unique.'
        )
      }
      
      # lista di variabili effettivamente disponibili in tutti gli mgnet (union)
      available <- meta_taxa_vars(object, .fmt = "unique")
      missing   <- setdiff(groups, available)
      if (length(missing)) {
        cli::cli_abort(
          c(
            "x" = 'Unknown grouping variable{?s} stored in attribute {.field "mgnet_groups"} on {.cls mgnets}: {.val {missing}}.',
            "v" = "Available variables across all elements are: {.val {available}}."
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
#' Build an `mgnets` object to encapsulate a collection of `mgnet` objects.
#' Accepts multiple `mgnet` arguments or a single *list* of `mgnet`.
#'
#' @param ... `mgnet` objects, or a *list* of `mgnet` objects.
#'
#' @return An object of class `mgnets`. If not empty, elements must be named with
#' unique names.
#'
#' @importFrom methods new
#' @export
#' @name mgnets
mgnets <- function(...) {
  x <- list(...)

  # single list argument
  if (length(x) == 1L && is.list(x[[1L]])) {
    x <- x[[1L]]
  }

  if (length(x) > 0L) {
    # type check
    if (!all(vapply(x, function(z) inherits(z, "mgnet"), logical(1)))) {
      cli::cli_abort("All inputs must be {.cls mgnet} objects or a list of {.cls mgnet}.")
    }
    # names present, non-empty
    nm <- names(x)
    if (is.null(nm) || any(!nzchar(nm))) {
      cli::cli_abort("All mgnet objects must be named with non-empty strings.")
    }
    # unique names
    if (anyDuplicated(nm) > 0L) {
      cli::cli_abort("Names must be unique.")
    }
  }

  methods::new("mgnets", mgnets = x)
}

#------------------------------------------------------------------------------#
# `mgnets` General methods 
#------------------------------------------------------------------------------#

#' General methods for `mgnets` objects
#'
#' Basic utilities to interact with `mgnets` objects: coerce to list, get length,
#' get and set names.
#'
#' @section Methods:
#' \describe{
#'   \item{`as.list(x)`}{Convert an `mgnets` object to a plain list of `mgnet` objects.}
#'   \item{`as(x, 'list')`}{Coerce the `mgnets` object to a list.}
#'   \item{`length(x)`}{Return the number of contained `mgnet` objects.}
#'   \item{`names(x)`}{Return element names.}
#'   \item{`names(x) <- value`}{Set element names.}
#' }
#'
#' @param x An object of class `mgnets`.
#' @param value A character vector for the replacement form `names<-`.
#' @param ... Not used.
#'
#' @return
#' - `as.list()` / `as(., 'list')`: a list of `mgnet` objects. \cr
#' - `length()`: an integer with the number of elements. \cr
#' - `names()`: a character vector with element names. \cr
#' - `names<-()`: the updated `mgnets` object.
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

# --- as.list (S3) ---
#' @rdname mgnets-methods
#' @export
as.list.mgnets <- function(x, ...) {
  x@mgnets
}

# --- Coercion to list (S4 setAs) ---
# no roxigen
methods::setAs("mgnets", "list", function(from) from@mgnets)

# --- length (S4) ---
#' @rdname mgnets-methods
#' @export
setMethod("length", "mgnets", function(x) {
  length(x@mgnets)
})

# --- names (getter S4) ---
#' @rdname mgnets-methods
#' @export
setMethod("names", "mgnets", function(x) {
  names(x@mgnets)
})

# --- names<- (replacement S4) ---
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
#' List-like accessors using `$`, `$<-`, `[[`, and `[[<-`.
#'
#' @param x An `mgnets` object.
#' @param name Length-1 character: the element name (for `$`).
#' @param i Index or name (for `[[` / `[[<-`).
#' @param j Ignored (signature compatibility).
#' @param value An `mgnet` object to assign.
#' @param ... Not used.
#'
#' @return
#' - `$` / `[[`: the selected `mgnet` object. \cr
#' - `$<-` / `[[<-`: the updated `mgnets` object.
#'
#' @aliases $,mgnets-method
#' $<-,mgnets,mgnet-method
#' [[,mgnets-method
#' [[<-,mgnets,ANY,ANY,mgnet-method
#'
#' @name mgnets-access
#' @rdname mgnets-access
NULL

# --- $ (getter S4) ---
#' @rdname mgnets-access
#' @export
setMethod("$", "mgnets", function(x, name) {
  if (!name %in% names(x@mgnets)) {
    cli::cli_abort("No {.cls mgnet} object named {.val {name}} found in the {.cls mgnets} object.")
  }
  x@mgnets[[name]]
})

# --- $<- (replacement S4) ---
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

# --- [[ (getter S4) ---
#' @rdname mgnets-access
#' @export
setMethod("[[", "mgnets", function(x, i, j, ...) {
  x@mgnets[[i]]
})

# --- [[<- (replacement S4) ---
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