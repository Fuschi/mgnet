# ABUNDANCE
#------------------------------------------------------------------------------#
#' Set Abundance Data
#'
#' This setter function allows you to update the abundance data for `mgnet` objects
#' and each `mgnet` object within a `mgnetList`. The abundance data must be a numeric matrix
#' for `mgnet` objects. For `mgnetList` objects, the abundance data should be a named list
#' of numeric matrices corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new abundance data to be set.
#' @return The `mgnet` or `mgnetList` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name abundance<-
#' @aliases abundance<-,mgnet-method abundance<-,mgnetList-method
setGeneric("abundance<-", function(object, value) standardGeneric("abundance<-"))

setMethod("abundance<-", c("mgnet","ANY"), function(object, value){
  object@abundance <- value
  validObject(object)
  object
})

setMethod("abundance<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@abundance <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# REL_ABUNDANCE
#------------------------------------------------------------------------------#
#' Set rel_abundance Data
#'
#' This setter function allows you to update the rel_abundance data for `mgnet` objects
#' and each `mgnet` object within a `mgnetList`. The rel_abundance data must be a numeric matrix
#' for `mgnet` objects. For `mgnetList` objects, the rel_abundance data should be a named list
#' of numeric matrices corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new rel_abundance data to be set.
#' @return The `mgnet` or `mgnetList` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name rel_abundance<-
#' @aliases rel_abundance<-,mgnet-method rel_abundance<-,mgnetList-method
setGeneric("rel_abundance<-", function(object, value) standardGeneric("rel_abundance<-"))

setMethod("rel_abundance<-", c("mgnet","ANY"), function(object, value){
  object@rel_abundance <- value
  validObject(object)
  object
})

setMethod("rel_abundance<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@rel_abundance <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# norm_abundance
#------------------------------------------------------------------------------#
#' Set Log-Transformed Abundance Data
#'
#' This setter function allows you to update the log-transformed abundance data for `mgnet` objects
#' and each `mgnet` object within a `mgnetList`. The log-transformed abundance data must be a numeric matrix
#' for `mgnet` objects. For `mgnetList` objects, the data should be a named list of numeric matrices
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new log-transformed abundance data to be set.
#' @return The `mgnet` or `mgnetList` object with updated log-transformed abundance data.
#' @export
#' @importFrom methods validObject
#' @name norm_abundance<-
#' @aliases norm_abundance<-,mgnet-method norm_abundance<-,mgnetList-method
setGeneric("norm_abundance<-", function(object, value) standardGeneric("norm_abundance<-"))

setMethod("norm_abundance<-", "mgnet", function(object, value) {
  object@norm_abundance <- value
  validObject(object)
  object
})

setMethod("norm_abundance<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@norm_abundance <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})



# INFO_SAMPLE
#------------------------------------------------------------------------------#
#' Set Sample Metadata
#'
#' This setter function updates the sample metadata for `mgnet` objects and
#' each `mgnet` object within a `mgnetList`. The metadata must be provided as a dataframe
#' for `mgnet` objects. For `mgnetList` objects, the metadata should be a named list of dataframes
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new sample metadata to be set.
#' @return The `mgnet` or `mgnetList` object with updated sample metadata.
#' @export
#' @importFrom methods validObject
#' @name info_sample<-
#' @aliases info_sample<-,mgnet-method info_sample<-,mgnetList-method
setGeneric("info_sample<-", function(object, value) standardGeneric("info_sample<-"))

setMethod("info_sample<-", c("mgnet","ANY"), function(object, value){
  object@info_sample <- value
  validObject(object)
  object
})

setMethod("info_sample<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@info_sample <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# LINEAGE
#------------------------------------------------------------------------------#
#' Set Lineage Data
#'
#' This setter function updates the lineage data for `mgnet` objects and
#' each `mgnet` object within a `mgnetList`. The lineage data must be provided as a character matrix
#' for `mgnet` objects. For `mgnetList` objects, the lineage data should be a named list of character matrices
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new lineage data to be set.
#' @return The `mgnet` or `mgnetList` object with updated lineage data.
#' @export
#' @importFrom methods validObject
#' @name lineage<-
#' @aliases lineage<-,mgnet-method lineage<-,mgnetList-method
setGeneric("lineage<-", function(object, value) standardGeneric("lineage<-"))

setMethod("lineage<-", "mgnet", function(object, value) {
  object@lineage <- value
  validObject(object)
  object
})

setMethod("lineage<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@lineage <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# INFO_TAXA
#------------------------------------------------------------------------------#
#' Set Taxa Metadata
#'
#' This setter function updates the taxa metadata for `mgnet` objects and
#' each `mgnet` object within a `mgnetList`. The taxa metadata must be provided as a dataframe
#' for `mgnet` objects. For `mgnetList` objects, the taxa metadata should be a named list of dataframes
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new taxa metadata to be set.
#' @return The `mgnet` or `mgnetList` object with updated taxa metadata.
#' @export
#' @importFrom methods validObject
#' @name info_taxa<-
#' @aliases info_taxa<-,mgnet-method info_taxa<-,mgnetList-method
setGeneric("info_taxa<-", function(object, value) standardGeneric("info_taxa<-"))

setMethod("info_taxa<-", "mgnet", function(object, value) {
  object@info_taxa <- value
  validObject(object)
  object
})

setMethod("info_taxa<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@info_taxa <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# NETWORK
#------------------------------------------------------------------------------#
#' Set Network Data
#'
#' This setter function updates the network data for `mgnet` objects and
#' each `mgnet` object within a `mgnetList`. The network data must be an `igraph` object
#' for `mgnet` objects. For `mgnetList` objects, the network data should be a named list of `igraph` objects
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new network data to be set.
#' @return The `mgnet` or `mgnetList` object with updated network data.
#' @export
#' @importFrom methods validObject
#' @name network<-
#' @aliases network<-,mgnet-method network<-,mgnetList-method
setGeneric("network<-", function(object, value) standardGeneric("network<-"))

setMethod("network<-", "mgnet", function(object, value) {
  object@network <- value
  validObject(object)
  object
})

setMethod("network<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@network <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})


# COMMUNITY
#------------------------------------------------------------------------------#
#' Set Community Detection Results
#'
#' This setter function updates the community detection results for `mgnet` objects and
#' each `mgnet` object within a `mgnetList`. The community data must be an object of class `communities`
#' for `mgnet` objects. For `mgnetList` objects, the community data should be a named list of `communities` objects
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new community detection results to be set, which should be compatible
#' with the structure expected by the `mgnet` object's community slot. For individual `mgnet`
#' objects, this is typically an object of class `communities` as returned by community detection
#' functions in the `igraph` package. For `mgnetList` objects, provide a named list where each element
#' is a `communities` object corresponding to the respective `mgnet` object within the list.
#'
#' @return The `mgnet` or `mgnetList` object with updated community detection results.
#'
#' @export
#' @importFrom methods validObject
#' @name community<-
#' @aliases community<-,mgnet-method community<-,mgnetList-method
setGeneric("community<-", function(object, value) standardGeneric("community<-"))

setMethod("community<-", "mgnet", function(object, value) {
  if(!inherits(value, "communities")) {
    stop("The value for 'community' must be a 'communities' object from the 'igraph' package.")
  }
  object@community <- value
  validObject(object)
  object
})

setMethod("community<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in 1:length(object)) {
    object@mgnets[[i]]@community <- value[[i]]
    validObject(object@mgnets[[i]])
  }
  validObject(object)
  object
})

