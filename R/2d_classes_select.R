# MODIFY AND UPDATE SAMPLE INFORMATION
#------------------------------------------------------------------------------#
#' Modify and Update Sample Information in mgnet or mgnetList Objects
#'
#' @description
#' Modifies the `meta` data frame within `mgnet` objects by selecting specific columns
#' using the `select` function from the `dplyr` package, and updates the `meta` slot with
#' the new data frame. This method provides dynamic manipulation of sample metadata,
#' allowing selection of columns based on specific analysis needs.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Columns to select from the `meta` data frame, specified unquoted due to
#'        non-standard evaluation supported by `dplyr`.
#' @return The modified `mgnet` or `mgnetList` object with an updated `meta` data frame.
#'
#' @details
#' Leverages the `select` function from `dplyr` for intuitive and flexible column selection,
#' ensuring that only relevant metadata are retained in the `meta` slot of the `mgnet` object.
#'
#' @export
#' @importFrom dplyr select
#' @name select_meta
#' @aliases select_meta,mgnet-method select_meta,mgnetList-method
setGeneric("select_meta", function(object, ...) standardGeneric("select_meta"))

setMethod("select_meta", signature(object = "mgnet"), function(object, ...) {
  new_meta <- dplyr::select(meta(object), ...)
  meta(object) <- new_meta
  object
})

setMethod("select_meta", signature(object = "mgnetList"), function(object, ...) {
  
  meta_new <- meta(object, .fmt = "tbl") %>%
    dplyr::select("mgnet", "sample_id", ...) %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  meta(object) <- meta_new
  object
})


# MODIFY AND UPDATE TAXA INFORMATION
#------------------------------------------------------------------------------#
#' Modify and Update taxa Information in mgnet or mgnetList Objects
#'
#' @description
#' Modifies the `taxa` data frame within `mgnet` objects by selecting specific columns
#' using the `select` function from the `dplyr` package, and updates the `taxa` slot with
#' the new data frame. This method provides dynamic manipulation of taxa taxadata,
#' allowing selection of columns based on specific analysis needs.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Columns to select from the `taxa` data frame, specified unquoted due to
#'        non-standard evaluation supported by `dplyr`.
#' @return The modified `mgnet` or `mgnetList` object with an updated `taxa` data frame.
#'
#' @details
#' Leverages the `select` function from `dplyr` for intuitive and flexible column selection,
#' ensuring that only relevant taxadata are retained in the `taxa` slot of the `mgnet` object.
#'
#' @export
#' @importFrom dplyr select
#' @name select_taxa
#' @aliases select_taxa,mgnet-method select_taxa,mgnetList-method
setGeneric("select_taxa", function(object, ...) standardGeneric("select_taxa"))

setMethod("select_taxa", signature(object = "mgnet"), function(object, ...) {
  new_taxa <- dplyr::select(taxa(object), ...)
  taxa(object) <- new_taxa
  object
})

setMethod("select_taxa", signature(object = "mgnetList"), function(object, ...) {
  
  taxa_new <- taxa(object, .fmt = "tbl") %>%
    dplyr::select("mgnet", "taxa_id", ...) %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  taxa(object) <- taxa_new
  object
})