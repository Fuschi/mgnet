#' #' Merge Two mgnet Objects by Samples
#' #'
#' #' @description
#' #' Merges two `mgnet` objects by aligning taxa, stacking their sample data vertically.
#' #'
#' #' @param mgnet1 First `mgnet` object.
#' #' @param mgnet2 Second `mgnet` object.
#' #' @return An `mgnet` object containing combined data from both input objects.
#' #' @export
#' setGeneric("bind_samples", function(mgnet1, mgnet2) {standardGeneric("bind_samples")})
#' 
#' setMethod("bind_samples", signature(c("mgnet", "mgnet")), function(mgnet1, mgnet2) {
#'   
#'   # CHECKS
#'   
#'   
#'   
#'   validObject(new_mgnet)
#'   return(new_mgnet)
#' })