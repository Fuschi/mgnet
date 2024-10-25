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
#' @importFrom purrr map imap
setGeneric("unite_meta", function(object, col, ..., sep = "_", remove = TRUE) {
  standardGeneric("unite_meta")
})

setMethod("unite_meta", "mgnet", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if meta data exists
  if(miss_sample(object)) {stop("Error: No sample available.")}
  if(length(meta(object)) == 0) stop("No meta data available in 'mgnet' object.")
  
  # Use dplyr to unite columns
  meta_unite <- tidyr::unite(meta(object, .fmt = "tbl"), {{col}}, ..., sep = sep, remove = FALSE)
  
  if(remove){
    cols_to_remove <- setdiff(sapply(rlang::ensyms(...), rlang::as_string), "sample_id")
    meta_unite <- dplyr::select(meta_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  meta(object) <- tibble::column_to_rownames(meta_unite, "sample_id")
  
  # Return the modified object
  return(object)
})

setMethod("unite_meta", "mgnetList", function(object, col, ..., sep = "_", remove = TRUE) {
  # Check if meta data exists
  if (any(length(meta(object)) == 0)) stop("No meta data available in 'mgnetList' object.")
  
  # Use dplyr to unite columns
  meta_merged <- tidyr::unite(meta(object, .fmt = "tbl"), {{col}}, ..., sep = sep, remove = FALSE) 
  
  if(remove){
    cols_to_remove <- setdiff(sapply(rlang::ensyms(...), rlang::as_string), c("mgnet","sample_id"))
    meta_merged <- dplyr::select(meta_merged, -tidyselect::any_of(cols_to_remove))
  }
  
  meta_splitted <- meta_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-tidyselect::all_of("mgnet")) %>%
        tibble::column_to_rownames("sample_id")
    })
  
  meta_splitted <- meta_splitted[names(object)]
  meta(object) <- meta_splitted
  
  # Return the modified object
  return(object)
})


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
  if(miss_taxa(object)) {stop("Error: No taxa available.")}
  if(length(taxa(object)) == 0) stop("No taxa data available in 'mgnet' object.")
  
  # Use dplyr to unite columns
  taxa_unite <- tidyr::unite(taxa(object, .fmt = "tbl"), {{col}}, ..., sep = sep, remove = FALSE) %>%
    dplyr::select(-tidyselect::any_of("comm_id"))
  
  if(remove){
    cols_to_remove <- setdiff(sapply(rlang::ensyms(...), rlang::as_string), "taxa_id")
    taxa_unite <- dplyr::select(taxa_unite, -tidyselect::any_of(cols_to_remove))
  }
  
  taxa(object) <- tibble::column_to_rownames(taxa_unite, "taxa_id")
  
  # Return the modified object
  return(object)
})

setMethod("unite_taxa", "mgnetList", function(object, col, ..., sep = "_", remove = TRUE) {
  if(miss_sample(object, "any")) {stop("Error: No sample available in any of the mgnet objects.")}
  if(miss_slot(object, "taxa", "all")) stop("No taxa data available in 'mgnetList' object.")
  
  # Use dplyr to unite columns
  taxa_merged <- tidyr::unite(taxa(object, .fmt = "tbl"), {{col}}, ..., sep = sep, remove = FALSE) 
  
  if(remove){
    cols_to_remove <- setdiff(sapply(rlang::ensyms(...), rlang::as_string), c("mgnet","taxa_id"))
    taxa_merged <- dplyr::select(taxa_merged, -tidyselect::any_of(cols_to_remove))
  }
  
  taxa_splitted <- taxa_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-tidyselect::any_of(c("mgnet", "comm_id"))) %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  taxa_splitted <- taxa_splitted[names(object)]
  taxa(object) <- taxa_splitted
  
  # Return the modified object
  return(object)
})