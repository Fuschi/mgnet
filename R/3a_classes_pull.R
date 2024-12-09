# PULL INFO SAMPLE
#------------------------------------------------------------------------------#
#' Pull Sample Information from mgnet or mgnetList Objects
#'
#' @description
#' Retrieves specific sample information from the `meta` data frame within `mgnet` objects,
#' or each `mgnet` object within an `mgnetList`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param var The name or position of the column to retrieve from the `meta` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `mgnet` object, a vector containing the data from the specified column 
#'         of the `meta` data frame is returned. For an `mgnetList` object, a list of such vectors
#'         is returned, each corresponding to one `mgnet` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names.
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom rlang ensym
#' @importFrom purrr map
#' @name pull_meta
#' @aliases pull_meta,mgnet-method pull_meta,mgnetList-method
setGeneric("pull_meta", function(object, var = -1) standardGeneric("pull_meta"))

setMethod("pull_meta", signature(object = "mgnet"), function(object, var = -1) {
  if(miss_sample(object)) {stop("Error: No sample available.")}
  var <- rlang::ensym(var)
  dplyr::pull(gather_meta(object), var)
})

setMethod("pull_meta", signature(object = "mgnetList"), function(object, var = -1) {
  if(miss_sample(object, "any")) {stop("Error: No sample available in any of the mgnet objects.")}
  var <- rlang::ensym(var)
  gather_meta(object) %>% dplyr::select(mgnet, var) %>% base::split(.[["mgnet"]]) %>%
    purrr::map(\(x) pull(x, !!var))
  
})


# PULL TAXA
#------------------------------------------------------------------------------#
#' Pull Taxa Information from mgnet or mgnetList Objects
#'
#' @description
#' Retrieves specific taxonomic information from the `taxa` data frame within `mgnet` objects,
#' or each `mgnet` object within an `mgnetList`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param var The name or position of the column to retrieve from the `taxa` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `mgnet` object, a vector containing the data from the specified column 
#'         of the `taxa` data frame is returned. For an `mgnetList` object, a list of such vectors
#'         is returned, each corresponding to one `mgnet` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names, allowing more
#' flexible and intuitive usage within data manipulation workflows.
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom rlang ensym
#' @importFrom purrr map
#' @name pull_taxa
#' @aliases pull_taxa,mgnet-method pull_taxa,mgnetList-method
setGeneric("pull_taxa", function(object, var = -1) standardGeneric("pull_taxa"))

setMethod("pull_taxa", signature(object = "mgnet"), function(object, var = -1) {
  if(miss_taxa(object)) {stop("Error: No taxa available.")}
  var <- rlang::ensym(var)
  dplyr::pull(gather_taxa(object), var)
})

setMethod("pull_taxa", signature(object = "mgnetList"), function(object, var = -1) {
  if(miss_taxa(object, "any")) {stop("Error: No taxa available in any of the mgnet objects.")}
  var <- rlang::ensym(var)
  gather_taxa(object) %>% dplyr::select(mgnet, var) %>% base::split(.[["mgnet"]]) %>%
    purrr::map(\(x) pull(x, !!var))
  
})