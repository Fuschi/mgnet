#' Unite Multiple Sample Metadata Columns into One in `mgnet` or `mgnetList` Objects
#'
#' @description
#' Combines multiple columns from the metadata (`meta`) of `mgnet` or `mgnetList` objects into a single column, 
#' creating a new column by concatenating the contents of specified columns separated by a specified separator.
#'
#' @param object An `mgnet` or `mgnetList` object containing the metadata in which columns are to be united.
#' @param col Character string specifying the name of the new column that will be created by uniting existing columns.
#' @param ... Columns to be united into the new column. These should be specified unquoted and separated by commas.
#' @param sep A string specifying the separator to be used to concatenate the values of the columns being united. 
#'        Default is an underscore ("_").
#' @param remove A logical value indicating whether the original columns should be removed after uniting. 
#'        Default is `TRUE`.
#'
#' @return Returns the `mgnet` or `mgnetList` object with updated metadata where specified columns are united into 
#'         a new single column. For `mgnetList`, the operation is performed for each `mgnet` object within the list.
#'
#' @details
#' The `unite_meta` function is a wrapper around `tidyr::unite` tailored for handling the metadata of `mgnet` objects. 
#' It allows for the consolidation of information spread across multiple columns into a single column, which can be useful 
#' for creating unique identifiers or combining categorical variables.
#'
#' For `mgnetList` objects, the function applies the unite operation to each `mgnet` object contained within the list, 
#' ensuring that all metadata across the list is consistently processed.
#'
#' @seealso 
#' \code{\link[tidyr::unite]{tidyr}} for more details on the underlying unite function.
#' \code{\link[mgnet]{mgnet}} and \code{\link[mgnet]{mgnetList}} for details on the data structures.
#'
#' @name unite_meta
#' @aliases unite_meta,mgnet-method unite_meta,mgnetList-method
#' @export
#' @importFrom tidyr unite
#' @importFrom tidyselect all_of
#' @importFrom dplyr select arrange
#' @importFrom tibble column_to_rownames
#' @importFrom purrr map imap map_chr
#' @importFrom rlang as_label enquos
setGeneric("unite_meta", function(object, col, ..., sep = "_", remove = TRUE) {
  standardGeneric("unite_meta")
})

setMethod("unite_meta", "mgnet", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if meta data exists
  if(miss_sample(object)) {stop("Error: No sample available.")}
  if(miss_slot(object, "meta")) stop("No meta data available in 'mgnet' object.")
  
  # Evaluate column selection to resolve tidyselect helpers
  meta_data <- gather_meta(object)
  selected_columns <- colnames(dplyr::select(meta_data, !!!rlang::enquos(...)))
  
  # Use tidyr::unite to combine columns
  meta_unite <- tidyr::unite(meta_data, !!col, !!!rlang::syms(selected_columns), sep = sep, remove = FALSE)

  if (remove) {
    cols_to_remove <- setdiff(selected_columns, "sample_id")
    meta_unite <- dplyr::select(meta_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  meta(object) <- meta_unite
  
  # Return the modified object
  return(object)
})

setMethod("unite_meta", "mgnetList", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if meta data exists
  if(miss_sample(object, "any")) {stop("Error: No sample available.")}
  if(miss_slot(object, "meta", "all")) stop("No meta data available in 'mgnet' object.")
  
  # Evaluate column selection to resolve tidyselect helpers
  meta_data <- gather_meta(object)
  selected_columns <- colnames(dplyr::select(meta_data, !!!rlang::enquos(...)))
  
  # Use tidyr::unite to combine columns
  meta_unite <- tidyr::unite(meta_data, !!col, !!!rlang::syms(selected_columns), sep = sep, remove = FALSE)
  
  if (remove) {
    cols_to_remove <- setdiff(selected_columns, c("mgnet","sample_id"))
    meta_unite <- dplyr::select(meta_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  meta(object) <- meta_unite
  
  # Return the modified object
  return(object)
})


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Unite Multiple Taxa Metadata Columns into One in `mgnet` or `mgnetList` Objects
#'
#' @description
#' Combines multiple columns from the metadata (`taxa`) of `mgnet` or `mgnetList` objects into a single column, 
#' creating a new column by concatenating the contents of specified columns separated by a specified separator.
#'
#' @param object An `mgnet` or `mgnetList` object containing the metadata in which columns are to be united.
#' @param col Character string specifying the name of the new column that will be created by uniting existing columns.
#' @param ... Columns to be united into the new column. These should be specified unquoted and separated by commas.
#' @param sep A string specifying the separator to be used to concatenate the values of the columns being united. 
#'        Default is an underscore ("_").
#' @param remove A logical value indicating whether the original columns should be removed after uniting. 
#'        Default is `TRUE`.
#'
#' @return Returns the `mgnet` or `mgnetList` object with updated metadata where specified columns are united into 
#'         a new single column. For `mgnetList`, the operation is performed for each `mgnet` object within the list.
#'
#' @details
#' The `unite_taxa` function is a wrapper around `tidyr::unite` tailored for handling the metadata of `mgnet` objects. 
#' It allows for the consolidation of information spread across multiple columns into a single column, which can be useful 
#' for creating unique identifiers or combining categorical variables.
#'
#' For `mgnetList` objects, the function applies the unite operation to each `mgnet` object contained within the list, 
#' ensuring that all metadata across the list is consistently processed.
#'
#' @seealso 
#' \code{\link[tidyr::unite]{tidyr}} for more details on the underlying unite function.
#' \code{\link[mgnet]{mgnet}} and \code{\link[mgnet]{mgnetList}} for details on the data structures.
#'
#' @name unite_taxa
#' @aliases unite_taxa,mgnet-method unite_taxa,mgnetList-method
#' @export
#' @importFrom tidyr unite
#' @importFrom tidyselect any_of all_of
#' @importFrom dplyr select arrange
#' @importFrom tibble column_to_rownames
#' @importFrom purrr map imap
setGeneric("unite_taxa", function(object, col, ..., sep = "_", remove = TRUE) {
  standardGeneric("unite_taxa")
})

setMethod("unite_taxa", "mgnet", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if taxa data exists
  if(miss_taxa(object)) {stop("Error: No sample available.")}
  if(miss_slot(object, "taxa")) stop("No taxa data available in 'mgnet' object.")
  
  # Evaluate column selection to resolve tidyselect helpers
  taxa_data <- gather_taxa(object)
  selected_columns <- colnames(dplyr::select(taxa_data, !!!rlang::enquos(...)))
  
  # Use tidyr::unite to combine columns
  taxa_unite <- tidyr::unite(taxa_data, !!col, !!!rlang::syms(selected_columns), sep = sep, remove = FALSE)
  
  if (remove) {
    cols_to_remove <- setdiff(selected_columns, "taxa_id")
    taxa_unite <- dplyr::select(taxa_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  taxa(object) <- taxa_unite
  
  # Return the modified object
  return(object)
})

setMethod("unite_taxa", "mgnetList", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if taxa data exists
  if(miss_sample(object, "any")) {stop("Error: No sample available.")}
  if(miss_slot(object, "taxa", "all")) stop("No taxa data available in 'mgnet' object.")
  
  # Evaluate column selection to resolve tidyselect helpers
  taxa_data <- gather_taxa(object)
  selected_columns <- colnames(dplyr::select(taxa_data, !!!rlang::enquos(...)))
  
  # Use tidyr::unite to combine columns
  taxa_unite <- tidyr::unite(taxa_data, !!col, !!!rlang::syms(selected_columns), sep = sep, remove = FALSE)
  
  if (remove) {
    cols_to_remove <- setdiff(selected_columns, c("mgnet","taxa_id"))
    taxa_unite <- dplyr::select(taxa_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  taxa(object) <- taxa_unite
  
  # Return the modified object
  return(object)
})
