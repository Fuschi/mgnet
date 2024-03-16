# ABUNDANCE
#------------------------------------------------------------------------------#
#' @title Set abundance for an `mgnet` object
#' 
#' @description Sets the abundance matrix in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A numeric matrix to set as the abundance matrix.
#' @usage abundance(object) <- value
#' @importFrom methods validObject
#' @export
"abundance<-" <- function(object, value) {
  ensure_mgnet(object)
  object@abundance <- value
  validObject(object)
  object
}

# LOG_ABUNDANCE
#------------------------------------------------------------------------------#
#' @title Set logabundance for an `mgnet` object
#' 
#' @description Sets the log-ratio transformed abundance matrix in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A numeric matrix to set as the log_abundance matrix.
#' @usage log_abundance(object) <- value
#' @importFrom methods validObject
#' @export
"log_abundance<-" <- function(object, value) {
  ensure_mgnet(object)
  object@log_abundance <- value
  validObject(object)
  object
}

# INFO_SAMPLE
#------------------------------------------------------------------------------#
#' @title Set info_sample for an `mgnet` object
#' 
#' @description Sets the info_sample matrix in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A data.frame to set as the sample metadata.
#' @usage info_sample(object) <- value
#' @importFrom methods validObject
#' @export
"info_sample<-" <- function(object, value) {
  ensure_mgnet(object)
  object@info_sample <- value
  validObject(object)
  object
}

# LINEAGE
#------------------------------------------------------------------------------#
#' @title Set lineage for an `mgnet` object
#' 
#' @description Sets the lineage matrix in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A character matrix to set as the taxonomic lineage matrix.
#' @usage lineage(object) <- value
#' @importFrom methods validObject
#' @export
"lineage<-" <- function(object, value) {
  ensure_mgnet(object)
  object@lineage <- value
  validObject(object)
  object
}

# INFO_TAXA
#------------------------------------------------------------------------------#
#' @title Set info_taxa for an `mgnet` object
#' 
#' @description Sets the info_taxa matrix in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A data.frame to set as the taxa metadata matrix.
#' @usage info_taxa(object) <- value
#' @importFrom methods validObject
#' @export
"info_taxa<-" <- function(object, value) {
  ensure_mgnet(object)
  object@info_taxa <- value
  validObject(object)
  object
}

# NETWORK
#------------------------------------------------------------------------------#
#' @title Set network for an `mgnet` object
#' 
#' @description Sets the `igraph` network in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value An `igraph` to set as the network.
#' @usage network(object) <- value
#' @importFrom methods validObject
#' @export
"network<-" <- function(object, value) {
  ensure_mgnet(object)
  object@network <- value
  validObject(object)
  object
}

# COMMUNITIES
#------------------------------------------------------------------------------#
#' @title Set communities for an `mgnet` object
#' 
#' @description Sets the communities in an `mgnet` object.
#'
#' @param object An object of class `mgnet`.
#' @param value A `communities` class object to set as the communities of the network.
#' @usage communities(object) <- value
#' @importFrom methods validObject
#' @export
"communities<-" <- function(object, value) {
  ensure_mgnet(object)
  object@communities <- value
  validObject(object)
  object
}