# REORDER SAMPLE
#------------------------------------------------------------------------------#
#' Reorder Samples in mgnet Object
#'
#' This function reorders the samples in an `mgnet` object according to a specified order.
#' The order can be specified using sample indices, sample identifiers, or a function that
#' returns the new order when applied to the object.
#'
#' @param object An `mgnet` object.
#' @param sample_order A numeric vector of indices, a character vector of sample labels, or
#'        a function that takes `object` as its input and returns a vector of new ordering indices.
#' @return The `mgnet` object with reordered samples.
#' @export
#' @name reorder_samples
#' @aliases reorder_samples,mgnet-method reorder_samples,mgnetList-method
setGeneric("reorder_samples", function(object, sample_order) standardGeneric("reorder_samples"))

setMethod("reorder_samples", "mgnet", function(object, sample_order) {
  
  # Evaluate if sample_order is a function and apply it
  if (is.function(sample_order)) {
    sample_order <- sample_order(object)
  }
  
  if (length(sample_order) != nsample(object) | !is.vector(sample_order)) {
    stop("sample_order must to be a vector with length equal to the number of samples.")
  }
  if (anyDuplicated(sample_order)) {
    stop("sample_order must contain unique entries.")
  }
  if(is.numeric(sample_order)){
    if(!all(1:nsample(object) %in% sample_order)) stop("the reordering must consider all the samples and someone missing")
  } else if(is.character(sample_order)){
    if(!all(sample_id(object) %in% sample_order)) stop("the reordering must consider all the samples and someone missing")
  } else {
    stop("sample order must to be a numeric or character vector that contain indices or sample names")
  }
  
  # Reorder
  object <- object[sample_order, ]
  
  # Additional reordering for other slots could be added here
  
  validObject(object)
  return(object)
})

setMethod("reorder_samples", "mgnetList", function(object, sample_order) {
  
  # Evaluate if sample_order is a function and apply it
  if (!is.function(sample_order)) {
    stop("sample_order for an mgnetList must to be a function only.")
  }
  
  for(i in 1:length(object)){
    indices_reordered_i <- sample_order(object@mgnets[[i]])
    object@mgnets[[i]] <- reorder_samples(object@mgnets[[i]], indices_reordered_i)
  }
  
  return(object)
})


# REORDER TAXA
#------------------------------------------------------------------------------#
#' Reorder Taxa in mgnet Object
#'
#' This function reorders the taxa in an `mgnet` object according to a specified order.
#' The order can be specified using taxa indices, taxa identifiers, or a function that
#' returns the new order when applied to the object.
#'
#' @param object An `mgnet` object.
#' @param taxa_order A numeric vector of indices, a character vector of taxa labels, or
#'        a function that takes `object` as its input and returns a vector of new ordering indices.
#' @return The `mgnet` object with reordered taxa.
#' @export
#' @name reorder_taxa
#' @aliases reorder_taxa,mgnet-method reorder_taxa,mgnetList-method
setGeneric("reorder_taxa", function(object, taxa_order) standardGeneric("reorder_taxa"))

setMethod("reorder_taxa", "mgnet", function(object, taxa_order) {
  
  # Evaluate if taxa_order is a function and apply it
  if (is.function(taxa_order)) {
    taxa_order <- taxa_order(object)
  }
  
  if (length(taxa_order) != ntaxa(object) | !is.vector(taxa_order)) {
    stop("taxa_order must to be a vector with length equal to the number of taxa.")
  }
  if (anyDuplicated(taxa_order)) {
    stop("taxa_order must contain unique entries.")
  }
  if(is.numeric(taxa_order)){
    if(!all(1:ntaxa(object) %in% taxa_order)) stop("the reordering must consider all the taxa and someone missing")
  } else if(is.character(taxa_order)){
    if(!all(taxa_id(object) %in% taxa_order)) stop("the reordering must consider all the taxa and someone missing")
  } else {
    stop("taxa order must to be a numeric or character vector that contain indices or taxa names")
  }
  
  # Reorder
  object <- object[, taxa_order]
  
  # Additional reordering for other slots could be added here
  
  validObject(object)
  return(object)
})

setMethod("reorder_taxa", "mgnetList", function(object, taxa_order) {
  
  # Evaluate if sample_order is a function and apply it
  if (!is.function(taxa_order)) {
    stop("taxa_order for an mgnetList must to be a function only.")
  }
  
  for(i in 1:length(object)){
    indices_reordered_i <- taxa_order(object@mgnets[[i]])
    object@mgnets[[i]] <- reorder_taxa(object@mgnets[[i]], indices_reordered_i)
  }
  
  return(object)
})