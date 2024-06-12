# SELECT INFO SAMPLE
#------------------------------------------------------------------------------#
#' Select Columns from info_sample in mgnet Objects
#'
#' This function allows users to select specific columns from the `info_sample` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's select semantics to provide a flexible interface for column selection based on column names or conditions.
#'
#' @description
#' `select_info_sample` uses dplyr's select functionality to enable precise selection of columns from the `info_sample`
#' data frame in `mgnet` objects. This can be particularly useful for simplifying metadata before further analysis
#' or for focusing on a subset of metadata attributes.
#'
#' @param object An `mgnet` or `mgnetList` object from which columns will be selected.
#' @param ... Conditions specifying which columns to select, passed to dplyr::select().
#'        This can include a variety of selectors like column names, dplyr helper functions (starts_with, ends_with, contains, etc.), or indices.
#'        For detailed usage, refer to \code{\link[dplyr]{select}}.
#'
#' @return An `mgnet` object with the `info_sample` slot updated to include only the selected columns, or an `mgnetList` object
#'         where each contained `mgnet` object has its `info_sample` slot similarly updated.
#'
#' @seealso
#' \code{\link[dplyr]{select}} for details on the selection conditions.
#' \code{\link{mgnet}} and \code{\link{mgnetList}} for details on the object structures.
#'
#' @export
#' @name select_info_sample
#' @aliases select_info_sample,mgnet-method select_info_sample,mgnetList-method
#' @importFrom dplyr select
#' @importFrom rlang enquos
setGeneric("select_info_sample", function(object, ...) {standardGeneric("select_info_sample")})

setMethod("select_info_sample", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  # Apply dplyr::select to info_sample
  selected_info_sample <- dplyr::select(info_sample(object, .fmt = "df"), !!!conditions)
  
  if(ncol(selected_info_sample) == 0) {
    stop("No columns meet the specified selecting conditions. Please revise your criteria.")
  }
  
  # Update info_sample directly
  object@info_sample <- selected_info_sample 
  validObject(object)
  
  return(object)
})


setMethod("select_info_sample", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    
    # Apply dplyr::select to info_sample
    selected_info_sample <- dplyr::select(info_sample(mgnet_obj, .fmt = "df"), !!!conditions)
    
    if(ncol(selected_info_sample) == 0) {
      stop("No columns meet the specified selecting conditions. Please revise your criteria.")
    }
    
    # Update info_sample directly
    mgnet_obj@info_sample <- selected_info_sample
    validObject(mgnet_obj)
    
    return(mgnet_obj)
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  validObject(object)
  return(object)
})

# SELECT INFO TAXA
#------------------------------------------------------------------------------#
#' Select Columns from info_taxa in mgnet Objects
#'
#' This function allows users to select specific columns from the `info_taxa` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's select semantics to provide a flexible interface for column selection based on column names or conditions.
#'
#' @description
#' `select_info_taxa` uses dplyr's select functionality to enable precise selection of columns from the `info_taxa`
#' data frame in `mgnet` objects. This can be particularly useful for simplifying metadata before further analysis
#' or for focusing on a subset of metadata attributes.
#'
#' @param object An `mgnet` or `mgnetList` object from which columns will be selected.
#' @param ... Conditions specifying which columns to select, passed to dplyr::select().
#'        This can include a variety of selectors like column names, dplyr helper functions (starts_with, ends_with, contains, etc.), or indices.
#'        For detailed usage, refer to \code{\link[dplyr]{select}}.
#'
#' @return An `mgnet` object with the `info_taxa` slot updated to include only the selected columns, or an `mgnetList` object
#'         where each contained `mgnet` object has its `info_taxa` slot similarly updated.
#'
#' @seealso
#' \code{\link[dplyr]{select}} for details on the selection conditions.
#' \code{\link{mgnet}} and \code{\link{mgnetList}} for details on the object structures.
#'
#' @export
#' @name select_info_taxa
#' @aliases select_info_taxa,mgnet-method select_info_taxa,mgnetList-method
#' @importFrom dplyr select
#' @importFrom rlang enquos
setGeneric("select_info_taxa", function(object, ...) {standardGeneric("select_info_taxa")})

setMethod("select_info_taxa", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  # Apply dplyr::select to info_taxa
  selected_info_taxa <- dplyr::select(info_taxa(object, .fmt = "df"), !!!conditions)
  
  if(ncol(selected_info_taxa) == 0) {
    stop("No columns meet the specified selecting conditions. Please revise your criteria.")
  }
  
  # Update info_sample directly
  object@info_taxa <- selected_info_taxa
  validObject(object)
  
  return(object)
})


setMethod("select_info_taxa", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    
    # Apply dplyr::select to info_sample
    selected_info_taxa <- dplyr::select(info_taxa(mgnet_obj, .fmt = "df"), !!!conditions)
    
    if(ncol(selected_info_taxa) == 0) {
      stop("No columns meet the specified selecting conditions. Please revise your criteria.")
    }
    
    # Update info_sample directly
    mgnet_obj@info_taxa <- selected_info_taxa
    validObject(mgnet_obj)
    
    return(mgnet_obj)
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  validObject(object)
  return(object)
})