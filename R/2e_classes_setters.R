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



# META
#------------------------------------------------------------------------------#
#' Update Sample Metadata for mgnet and mgnetList Objects
#'
#' This function sets the sample metadata for `mgnet` objects and also updates
#' each `mgnet` object within a `mgnetList`. The metadata for `mgnet` objects must be
#' provided as a dataframe. For `mgnetList` objects, the metadata should be either
#' a named list of dataframes corresponding to each `mgnet` object within the list,
#' or a single dataframe containing a 'mgnet' and 'sample_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `mgnet` or `mgnetList` object to be updated.
#' @param value The new sample metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `mgnet` or `mgnetList` object with new sample metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames has_rownames
#' @name meta<-
#' @aliases meta<-,mgnet-method meta<-,mgnetList-method
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))

setMethod("meta<-", c("mgnet", "ANY"), function(object, value) {
  
  if ("sample_id" %in% colnames(value) && tibble::has_rownames(value)) {
    stop(sprintf("Error: %s cannot have both the 'sample_id' column and row names set.", deparse(substitute(value))))
  }

  if (!"sample_id" %in% colnames(value) && !has_rownames(value)) {
    stop(sprintf("Error: %s must have either a 'sample_id' column or row names set.", deparse(substitute(value))))
  }

  if ("sample_id" %in% colnames(value)) {
    if (any(duplicated(value$sample_id))) {
      stop(sprintf("Error: The `sample_id` column in %s contains duplicate values.", deparse(substitute(value))))
    }
    value <- tibble::column_to_rownames(value, "sample_id")
  }

  if (has_sample(object)) {
    if (!all(rownames(value) %in% sample_id(object))) {
      stop(sprintf("Error: In '%s', the 'sample_id' values in the provided data do not match the existing 'sample_id' in the object.", i))
    }
    if (!identical(rownames(value), sample_id(object))) {
      value <- value[sample_id(object), ]
    }
  }

  object@meta <- as.data.frame(value)
  validObject(object)
  object
})

setMethod("meta<-", c("mgnetList","ANY"), function(object, value){
  
  if(class(value)[[1]] == "list"){
    
    are_list_assign(object, value)
    for(i in names(object)) meta(object[[i]]) <- value[[i]] 
    
  } else if(is.data.frame(value)){
    
    is_assign_tbl(object, value, "sample")
    splitted_value <- split(value, value$mgnet)
    splitted_value <- lapply(splitted_value, \(x){
      x$mgnet <- NULL
      return(x)
    })
    for(i in names(object)) meta(object[[i]]) <- splitted_value[[i]]
    
  } else {
    
    valueName <- deparse(substitute(value))
    stop(sprintf("Error: %s must could be only a list of data.frame or a single data.frame but with `mgnet` and `sample_id` columns"))
    
  }
  
  validObject(object)
  object
})


# TAXA
#------------------------------------------------------------------------------#
#' Update Taxa Metadata for mgnet and mgnetList Objects
#'
#' This function sets the taxa metadata for `mgnet` objects and also updates
#' each `mgnet` object within a `mgnetList`. The metadata for `mgnet` objects must be
#' provided as a dataframe. For `mgnetList` objects, the metadata should be either
#' a named list of dataframes corresponding to each `mgnet` object within the list,
#' or a single dataframe containing a 'mgnet' and 'taxa_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `mgnet` or `mgnetList` object to be updated.
#' @param value The new taxa metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `mgnet` or `mgnetList` object with new taxa metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames
#' @name taxa<-
#' @aliases taxa<-,mgnet-method taxa<-,mgnetList-method
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))

setMethod("taxa<-", c("mgnet", "ANY"), function(object, value) {
  
  if ("taxa_id" %in% colnames(value) && has_rownames(value)) {
    stop(sprintf("Error: %s cannot have both the 'taxa_id' column and row names set.", deparse(substitute(value))))
  }
  
  if (!"taxa_id" %in% colnames(value) && !has_rownames(value)) {
    stop(sprintf("Error: %s must have either a 'taxa_id' column or row names set.", deparse(substitute(value))))
  }
  
  if ("taxa_id" %in% colnames(value)) {
    if (any(duplicated(value$taxa_id))) {
      stop(sprintf("Error: The `taxa_id` column in %s contains duplicate values.", deparse(substitute(value))))
    }
    value <- tibble::column_to_rownames(value, "taxa_id")
  }
  
  if (has_sample(object)) {
    if (!all(rownames(value) %in% taxa_id(object))) {
      stop(sprintf("Error: In '%s', the 'taxa_id' values in the provided data do not match the existing 'taxa_id' in the object.", i))
    }
    if (!identical(rownames(value), taxa_id(object))) {
      value <- value[taxa_id(object), ]
    }
  }
  
  if ("comm_id" %in% colnames(value)){
    if(!all(value$comm_id == comm_id(object))){
      stop("Error: Communities membership encoded in comm_id differes from the ones in the object")
    }
    value$comm_id <- NULL
  }
  
  object@taxa <- as.data.frame(value)
  validObject(object)
  object
})

setMethod("taxa<-", c("mgnetList","ANY"), function(object, value){
  
  if(class(value)[[1]] == "list"){
    
    are_list_assign(object, value)
    for(i in names(object)) taxa(object[[i]]) <- value[[i]] 
    
  } else if(is.data.frame(value)){
    
    is_assign_tbl(object, value, "taxa")
    splitted_value <- split(value, value$mgnet)
    splitted_value <- lapply(splitted_value, \(x){
      x$mgnet <- NULL
      return(x)
    })
    for(i in names(object)) taxa(object[[i]]) <- splitted_value[[i]]
    
  } else {
    
    valueName <- deparse(substitute(value))
    stop(sprintf("Error: %s must could be only a list of data.frame or a single data.frame but with `mgnet` and `taxa_id` columns"))
    
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


# LINK
#------------------------------------------------------------------------------#
setGeneric("link<-", function(object, value, .suffix = c("_1", "_2")) standardGeneric("link<-"))

setMethod("link<-", c("mgnet", "ANY"), function(object, value) {
  
  if(miss_slot(object, "netw")) stop("Error: No network available.")
  
  correct_taxa_ids <- paste0("taxa", .suffix)
  col_names <- colnames(value)
  if(all(correct_taxa_ids %in% col_names)) stop("I didn't find the identifiers of the taxa for the nodes....")
  
  tbl_link <- value %>%
    dplyr::select(-tidyselect::ends_with(.suffix[1]),
                  -tidyselect::ends_with(.suffix[2]))
  
})

setMethod("taxa<-", c("mgnetList","ANY"), function(object, value){
  
  if(class(value)[[1]] == "list"){
    
    are_list_assign(object, value)
    for(i in names(object)) taxa(object[[i]]) <- value[[i]] 
    
  } else if(is.data.frame(value)){
    
    is_assign_tbl(object, value, "taxa")
    splitted_value <- split(value, value$mgnet)
    splitted_value <- lapply(splitted_value, \(x){
      x$mgnet <- NULL
      return(x)
    })
    for(i in names(object)) taxa(object[[i]]) <- splitted_value[[i]]
    
  } else {
    
    valueName <- deparse(substitute(value))
    stop(sprintf("Error: %s must could be only a list of data.frame or a single data.frame but with `mgnet` and `taxa_id` columns"))
    
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

