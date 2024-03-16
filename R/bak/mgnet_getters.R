# ABUNDANCE
#------------------------------------------------------------------------------#
#' @title Get abundance from an `mgnet`
#' 
#' @description
#' Retrieves the abundance matrix stored in an `mgnet` object. 
#'
#' @param object An object of class `mgnet`.
#' @return A numeric matrix representing the abundance matrix.
#' @usage abundance(object)
#' @export
abundance  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@abundance
}

# LOG_ABUNDANCE
#------------------------------------------------------------------------------#
#' @title Get log_abundance from an `mgnet`
#' 
#' @description
#' Retrieves the log-ratio transformed abundance matrix stored in an `mgnet` object. 
#'
#' @param object An object of class `mgnet`.
#' @return A numeric matrix representing the log_abundance matrix.
#' @usage log_abundance(object)
#' @export
log_abundance  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@log_abundance
}

# INFO_SAMPLE
#------------------------------------------------------------------------------#
#' @title Get info_sample from an `mgnet`
#' 
#' @description
#' Retrieves the sample metadata stored in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @return A data.frame representing the sample metadata.
#' @usage info_sample(object)
#' @export
info_sample  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@info_sample
}

# LINEAGE
#------------------------------------------------------------------------------#
#' @title Get lineage from an `mgnet`
#' 
#' @description
#' Retrieves the taxonomic lineage matrix stored in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @return A character matrix representing the taxonomic lineage.
#' @usage lineage(object)
#' @export
lineage  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@lineage
}

# INFO_TAXA
#------------------------------------------------------------------------------#
#' @title Get info_taxa from an `mgnet`
#' 
#' @description
#' Retrieves the taxa metadata stored in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @return A data.frame representing the taxa metadata.
#' @usage info_taxa(object)
#' @export
info_taxa  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@info_taxa
}

# NETWORK
#------------------------------------------------------------------------------#
#' @title Get network from an `mgnet`
#' 
#' @description
#' Retrieves the network (as an igraph object) stored in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @return An igraph object representing the taxa interaction network.
#' @usage network(object)
#' @export
network  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@network
}

# COMMUNITIES
#------------------------------------------------------------------------------#
#' @title Get communities from an `mgnet`
#' 
#' @description
#' Retrieves the communities detected within the network, stored in an `mgnet` object. 
#'
#' @param object An object of class `mgnet`.
#' @return An `communities` class object representing the community detection results.
#' @usage communities(object)
#' @export
communities  <- function(object) {
  if(!inherits(object, "mgnet")) {stop(sprintf("Input '%s' is not of class 'mgnet'.", deparse(substitute(object))))}
  object@communities
}