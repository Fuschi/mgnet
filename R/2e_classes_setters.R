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
#' @name abun<-
#' @aliases abun<-,mgnet-method abun<-,mgnetList-method
setGeneric("abun<-", function(object, value) standardGeneric("abun<-"))

setMethod("abun<-", c("mgnet","ANY"), function(object, value){
  object@abun <- value
  validObject(object)
  object
})

setMethod("abun<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@abun <- value[[i]] }
  validObject(object)
  object
})


# rela
#------------------------------------------------------------------------------#
#' Set rela Data
#'
#' This setter function allows you to update the rela data for `mgnet` objects
#' and each `mgnet` object within a `mgnetList`. The rela data must be a numeric matrix
#' for `mgnet` objects. For `mgnetList` objects, the rela data should be a named list
#' of numeric matrices corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new rela data to be set.
#' @return The `mgnet` or `mgnetList` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name rela<-
#' @aliases rela<-,mgnet-method rela<-,mgnetList-method
setGeneric("rela<-", function(object, value) standardGeneric("rela<-"))

setMethod("rela<-", c("mgnet","ANY"), function(object, value){
  object@rela <- value
  validObject(object)
  object
})

setMethod("rela<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@rela <- value[[i]] }
  validObject(object)
  object
})


# norm
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
#' @name norm<-
#' @aliases norm<-,mgnet-method norm<-,mgnetList-method
setGeneric("norm<-", function(object, value) standardGeneric("norm<-"))

setMethod("norm<-", "mgnet", function(object, value) {
  object@norm <- value
  validObject(object)
  object
})

setMethod("norm<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@norm <- value[[i]] }
  validObject(object)
  object
})



# meta
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
#' @name meta<-
#' @aliases meta<-,mgnet-method meta<-,mgnetList-method
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))

setMethod("meta<-", c("mgnet","ANY"), function(object, value){
  object@meta <- value
  validObject(object)
  object
})

setMethod("meta<-", c("mgnetList","ANY"), function(object, value){
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@meta <- value[[i]] }
  validObject(object)
  object
})


# taxa
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
#' @name taxa<-
#' @aliases taxa<-,mgnet-method taxa<-,mgnetList-method
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))

setMethod("taxa<-", "mgnet", function(object, value) {
  object@taxa <- value
  validObject(object)
  object
})

setMethod("taxa<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@taxa <- value[[i]] }
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
#' @name netw<-
#' @aliases netw<-,mgnet-method netw<-,mgnetList-method
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))

setMethod("netw<-", "mgnet", function(object, value) {
  object@netw <- value
  validObject(object)
  object
})

setMethod("netw<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@netw <- value[[i]] }
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
#' @name comm<-
#' @aliases comm<-,mgnet-method comm<-,mgnetList-method
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))

setMethod("comm<-", "mgnet", function(object, value) {
  object@comm <- value
  validObject(object)
  object
})

setMethod("comm<-", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@comm <- value[[i]] }
  validObject(object)
  object
})

