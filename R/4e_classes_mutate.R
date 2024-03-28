# MUTATE INFO SAMPLE
#------------------------------------------------------------------------------#
#' Mutate Columns in info_sample of mgnet Objects
#'
#' This function allows users to modify existing columns or add new columns to the `info_sample` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's mutate semantics for flexible data transformation.
#'
#' @description
#' `mutate_info_sample` leverages dplyr's mutate functionality to enable column transformations or additions within the `info_sample`
#' data frame in `mgnet` objects.
#'
#' @param object An `mgnet` or `mgnetList` object to be modified.
#' @param ... Transformations to apply, passed to dplyr::mutate().
#'
#' @return The `mgnet` or `mgnetList` object with its `info_sample` slot updated.
#'
#' @seealso
#' \code{\link[dplyr]{mutate}} for details on transformation conditions.
#'
#' @export
#' @aliases mutate_info_sample,mgnet-method mutate_info_sample,mgnetList-method
#' @importFrom dplyr mutate
setGeneric("mutate_info_sample", function(object, ...) {standardGeneric("mutate_info_sample")})

setMethod("mutate_info_sample", "mgnet", function(object, ...) {
  conditions <- rlang::enquos(...)
  mutated_info_sample <- dplyr::mutate(info_sample(object, .fmt = "df"), !!!conditions)
  object@info_sample <- mutated_info_sample 
  validObject(object)
  return(object)
})

setMethod("mutate_info_sample", "mgnetList", function(object, ...) {
  conditions <- rlang::enquos(...)
  object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
    mutated_info_sample <- dplyr::mutate(info_sample(mgnet_obj, .fmt = "df"), !!!conditions)
    mgnet_obj@info_sample <- mutated_info_sample
    validObject(mgnet_obj)
    return(mgnet_obj)
  })
  validObject(object)
  return(object)
})


# MUTATE INFO TAXA
#------------------------------------------------------------------------------#
#' Mutate Columns in info_taxa of mgnet Objects
#'
#' This function allows users to modify existing columns or add new columns to the `info_taxa` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's mutate semantics for flexible data transformation.
#'
#' @description
#' `mutate_info_taxa` leverages dplyr's mutate functionality to enable column transformations or additions within the `info_taxa`
#' data frame in `mgnet` objects. This can be particularly useful for adding new taxonomic information or adjusting existing ones.
#'
#' @param object An `mgnet` or `mgnetList` object to be modified.
#' @param ... Transformations to apply, passed to dplyr::mutate().
#' @return The `mgnet` or `mgnetList` object with its `info_taxa` slot updated.
#' @seealso
#' \code{\link[dplyr]{mutate}} for details on transformation conditions.
#' @export
#' @aliases mutate_info_taxa,mgnet-method mutate_info_taxa,mgnetList-method
#' @importFrom dplyr mutate
setGeneric("mutate_info_taxa", function(object, ...) {standardGeneric("mutate_info_taxa")})

setMethod("mutate_info_taxa", "mgnet", function(object, ...) {
  conditions <- rlang::enquos(...)
  mutated_info_taxa <- dplyr::mutate(info_taxa(object, .fmt = "df"), !!!conditions)
  object@info_taxa <- mutated_info_taxa 
  validObject(object)
  return(object)
})

setMethod("mutate_info_taxa", "mgnetList", function(object, ...) {
  conditions <- rlang::enquos(...)
  object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
    mutated_info_taxa <- dplyr::mutate(info_taxa(mgnet_obj, .fmt = "df"), !!!conditions)
    mgnet_obj@info_taxa <- mutated_info_taxa
    validObject(mgnet_obj)
    return(mgnet_obj)
  })
  validObject(object)
  return(object)
})
