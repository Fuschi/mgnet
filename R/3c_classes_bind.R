# MODIFY AND UPDATE SAMPLE INFORMATION BY BINDING COLUMNS
#------------------------------------------------------------------------------#
#' Bind New Columns to Sample Information in mgnet or mgnetList Objects
#'
#' @description
#' Adds new columns to the `meta` data frame within `mgnet` objects by binding them
#' using the `bind_cols` function from the `dplyr` package. It updates the `meta` slot with
#' the new data frame. This method allows you to augment your sample metadata with additional
#' variables, enhancing the information contained within your `mgnet` object.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Data frames, tibbles, or vectors to bind as new columns to the `meta` data frame.
#'        For `mgnet` objects, the additional columns should match the number of rows in the `meta` data frame.
#'        For `mgnetList` objects, the additional columns should be provided as named lists,
#'        where each name corresponds to an `mgnet` object in the list, and each element is a data frame,
#'        tibble, or vector to bind to that `mgnet` object's `meta` data frame.
#' @return The modified `mgnet` or `mgnetList` object with an updated `meta` data frame.
#'
#' @details
#' Leverages the `bind_cols` function from `dplyr` for intuitive and flexible column binding,
#' allowing you to augment the sample metadata with additional data.
#'
#' @importFrom dplyr bind_cols
#' @importFrom tidyselect all_of
#' @export
#' @name bind_meta
#' @aliases bind_meta,mgnet-method bind_meta,mgnetList-method
setGeneric("bind_meta", function(object, ...) standardGeneric("bind_meta"))

setMethod("bind_meta", signature(object = "mgnet"), function(object, ...) {
  if(miss_sample(object)) {stop("Error: No sample available.")}
  new_meta <- dplyr::bind_cols(object@meta, ...)
  meta(object) <- new_meta
  object
})

setMethod("bind_meta", signature(object = "mgnetList"), function(object, ...) {
  
  if(miss_sample(object, "any")) {stop("Error: No sample available in any of the mgnet objects.")}
  
  additional_args <- list(...)
  
  # Check if additional columns are provided as a named list matching mgnet names
  if (length(additional_args) == 1 && is.list(additional_args[[1]]) && all(names(additional_args[[1]]) %in% names(object@mgnets))) {
    # Additional columns are provided as a named list
    additional_cols_list <- additional_args[[1]]
    
    # Apply to each mgnet in the mgnetList
    object@mgnets <- purrr::imap(object@mgnets, function(mgnet_obj, mgnet_name) {
      if (mgnet_name %in% names(additional_cols_list)) {
        new_meta <- dplyr::bind_cols(meta(mgnet_obj), additional_cols_list[[mgnet_name]])
        meta(mgnet_obj) <- new_meta
      }
      mgnet_obj
    })
    
  } else {
    # Additional columns are to be bound to all mgnet objects' meta data
    # Get combined meta data
    meta_tbl <- meta(object, .fmt = "tbl")
    
    # Bind additional columns to the combined meta data
    new_meta_tbl <- dplyr::bind_cols(meta_tbl, ...)
    
    # Split back into individual mgnet meta data
    meta_new <- new_meta_tbl %>%
      split(.$mgnet) %>%
      purrr::imap(function(df, mgnet_name) {
        df %>%
          dplyr::arrange(match(sample_id, sample_id(object@mgnets[[mgnet_name]]))) %>%
          dplyr::select(-tidyselect::all_of(mgnet)) %>%
          tibble::column_to_rownames("sample_id")
      })
    
    # Assign the new meta data back to each mgnet
    meta(object) <- meta_new
  }
  
  object
})


# MODIFY AND UPDATE TAXA INFORMATION BY BINDING COLUMNS
#------------------------------------------------------------------------------#
#' Bind New Columns to Taxa Information in mgnet or mgnetList Objects
#'
#' @description
#' Adds new columns to the `taxa` data frame within `mgnet` objects by binding them
#' using the `bind_cols` function from the `dplyr` package. It updates the `taxa` slot with
#' the new data frame. This method allows you to augment your taxa metadata with additional
#' variables, enhancing the information contained within your `mgnet` object.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Data frames, tibbles, or vectors to bind as new columns to the `taxa` data frame.
#'        For `mgnet` objects, the additional columns should match the number of rows in the `taxa` data frame.
#'        For `mgnetList` objects, the additional columns should be provided as named lists,
#'        where each name corresponds to an `mgnet` object in the list, and each element is a data frame,
#'        tibble, or vector to bind to that `mgnet` object's `taxa` data frame.
#' @return The modified `mgnet` or `mgnetList` object with an updated `taxa` data frame.
#'
#' @details
#' Leverages the `bind_cols` function from `dplyr` for intuitive and flexible column binding,
#' allowing you to augment the taxa metadata with additional data.
#'
#' @importFrom dplyr bind_cols
#' @export
#' @name bind_taxa
#' @aliases bind_taxa,mgnet-method bind_taxa,mgnetList-method
setGeneric("bind_taxa", function(object, ...) standardGeneric("bind_taxa"))

setMethod("bind_taxa", signature(object = "mgnet"), function(object, ...) {
  if(miss_taxa(object)) {stop("Error: No taxa available.")}
  new_taxa <-  dplyr::bind_cols(object@taxa, ...)
  taxa(object) <- new_taxa
  object
})

setMethod("bind_taxa", signature(object = "mgnetList"), function(object, ...) {
  
  if(miss_taxa(object, "any")) {stop("Error: No taxa available in any of the mgnet objects.")}
  
  additional_args <- list(...)
  
  # Check if additional columns are provided as a named list matching mgnet names
  if (length(additional_args) == 1 && is.list(additional_args[[1]]) && all(names(additional_args[[1]]) %in% names(object@mgnets))) {
    # Additional columns are provided as a named list
    additional_cols_list <- additional_args[[1]]
    
    # Apply to each mgnet in the mgnetList
    object@mgnets <- purrr::imap(object@mgnets, function(mgnet_obj, mgnet_name) {
      if (mgnet_name %in% names(additional_cols_list)) {
        new_taxa <- dplyr::bind_cols(taxa(mgnet_obj), additional_cols_list[[mgnet_name]])
        taxa(mgnet_obj) <- new_taxa
      }
      mgnet_obj
    })
    
  } else {
    # Additional columns are to be bound to all mgnet objects' taxa data
    # Get combined taxa data
    taxa_tbl <- taxa(object, .fmt = "tbl")
    
    # Bind additional columns to the combined taxa data
    new_taxa_tbl <- dplyr::bind_cols(taxa_tbl, ...)
    
    # Split back into individual mgnet taxa data
    taxa_new <- new_taxa_tbl %>%
      split(.$mgnet) %>%
      purrr::imap(function(df, mgnet_name) {
        df %>%
          dplyr::arrange(match(taxa_id, taxa_id(object@mgnets[[mgnet_name]]))) %>%
          dplyr::select(-tidyselect::any_of(mgnet, comm_id)) %>%
          tibble::column_to_rownames("taxa_id")
      })
    
    # Assign the new taxa data back to each mgnet
    taxa(object) <- taxa_new
  }
  
  object
})
