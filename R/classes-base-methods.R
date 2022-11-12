setOldClass("ngs_data")
setOldClass("sample_metadata")
setOldClass("taxonomy_table")
setOldClass("mg")
################################################################################
################################################################################
# NSAMPLE
################################################################################
################################################################################
#' Get number of samples.
#' 
#' @description 
#' Return an integer indicating the number of sample.
#'
#' @usage data(object)
#'
#' @param object (Required) \code{\link{mg-class}} or \code{\link{ngs_data-class}}
#' or \code{\link{ngs_data-class}}.
#'
#' @rdname nsample-methods
#' @docType methods
#' @export
#' @aliases nsample nsample
setGeneric("nsample", function(object) standardGeneric("nsample"))
#' @rdname nsample-methods
#' @aliases nsample,ANY-method
setMethod("nsample", "ANY", function(object){NULL})
#' @rdname nsample-methods
#' @aliases nsample,ngs_data-method
setMethod("nsample", "ngs_data", function(object){
  return(ifelse(is.null(object),NULL,nrow(object@value)))
  })
#' @rdname nsample-methods
#' @aliases nsample,sample_metadata-method
setMethod("nsample", "sample_metadata", function(object){
  return(ifelse(is.null(object),NULL,nrow(object@value)))
})
#' @rdname nsample-methods
#' @aliases nsample,mg-method
setMethod("nsample", "mg", function(object){
  max(nsample(object@data),nsample(object@meta))
})
################################################################################
################################################################################
# NTAXA
################################################################################
################################################################################
#' Get number of taxa
#' 
#' @description 
#' Return an integer indicating the number of taxa.
#'
#' @usage taxa(object)
#'
#' @param object (Required) \code{\link{mg-class}} or \code{\link{ngs_data-class}}
#' or \code{\link{ngs_data-class}}.
#'
#' @rdname nsample-methods
#' @docType methods
#' @export
#' @aliases nsample nsample
setGeneric("nsample", function(object) standardGeneric("nsample"))
#' @rdname nsample-methods
#' @aliases nsample,ANY-method
setMethod("nsample", "ANY", function(object){NULL})
#' @rdname nsample-methods
#' @aliases nsample,ngs_data-method
setMethod("nsample", "ngs_data", function(object){
  return(ifelse(is.null(object),NULL,nrow(object@value)))
})
#' @rdname nsample-methods
#' @aliases nsample,sample_metadata-method
setMethod("nsample", "sample_metadata", function(object){
  return(ifelse(is.null(object),NULL,nrow(object@value)))
})
#' @rdname nsample-methods
#' @aliases nsample,mg-method
setMethod("nsample", "mg", function(object){
  max(nsample(object@data),nsample(object@meta))
})