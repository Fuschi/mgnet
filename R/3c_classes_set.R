# SET ABUNDANCE
#------------------------------------------------------------------------------#
#' Set `abundance` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the abundance data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new abundance data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated abundance data.
#' @export
#' @name set_abundance
#' @aliases set_abundance,mgnet-method set_abundance,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_abundance", function(object, value) standardGeneric("set_abundance"))

setMethod("set_abundance", "mgnet", function(object, value) {
  object@abundance <- value
  validObject(object)
  return(object)
})

setMethod("set_abundance", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_abundance(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})


# SET LOG_ABUNDANCE
#------------------------------------------------------------------------------#
#' Set `log_abundance` Slot in `mgnet` Objects
#'
#' This function allows setting the log_abundance data directly or calculating it
#' from the abundance data using specified methods ('clr' or 'iclr') after applying
#' zero replacement strategies ('unif' or 'const'). The zero replacement is applied
#' to the abundance data before the log-ratio transformation when 'method' is specified.
#'
#' @param object An `mgnet` object.
#' @param value Optional. A numeric matrix to set as log_abundance data directly.
#' @param clr_variant Optional. A string specifying the method to calculate log_abundance ('clr' or 'iclr').
#' @param zero_strategy A string specifying the zero replacement strategy applied before log-ratio transformation.
#' @return The `mgnet` object with updated log_abundance data.
#' @export
#' @seealso \code{\link{clr}}, \code{\link{iclr}}, \code{\link{zero_dealing}}
#' @name set_log_abundance
#' @aliases set_log_abundance,mgnet-method set_log_abundance,mgnetList-method
setGeneric("set_log_abundance", function(object,
                                         value = NULL,
                                         clr_variant = NULL, zero_strategy = "unif") standardGeneric("set_log_abundance"))

setMethod("set_log_abundance", "mgnet", function(object, value = NULL, clr_variant = NULL, zero_strategy = "unif") {

  if(is.null(value) && is.null(clr_variant)) {
    stop("Either 'value' or 'clr_variant' must be provided, not both.")
  }
  
  if(is.null(value) && is.null(clr_variant)) {
    stop("Either 'value' or 'clr_variant' must be provided.")
  }
  
  if(!is.null(value)) {
    
    object@log_abundance <- value
    validObject(object)
    return(object)
    
  }
  
  if(!is.null(clr_variant) && !is.character(clr_variant)) stop("`clr_variant` if setted it must be a string with possible choices 'clr' and 'iclr'")
  clr_variant <- match.arg(clr_variant, c("clr","iclr"))
  zero_strategy <- match.arg(zero_strategy, c("unif", "const"))
  if(length(object@abundance) == 0) stop("Abundance matrix missing; cannot calculate log-ratio abundance.")
  
  abundance_nozero <- zero_dealing(object@abundance, method = zero_strategy)
  
  if(clr_variant == "clr") {
    object@log_abundance <- clr(abundance_nozero)
  } else {
    object@log_abundance <- iclr(abundance_nozero)
  }
  
  validObject(object)
  return(object)
  
})


setMethod("set_log_abundance", "mgnetList", function(object, value = NULL,
                                                     clr_variant = NULL, zero_strategy = "unif") {
  
  if(!is.null(value)) are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_log_abundance(object@mgnets[[i]], value[[i]], clr_variant, zero_strategy)
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
  
})

# SET INFO_SAMPLE
#------------------------------------------------------------------------------#
#' Set `info_sample` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the info_sample data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new info_sample data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated info_sample data.
#' @export
#' @name set_info_sample
#' @aliases set_info_sample,mgnet-method set_info_sample,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_info_sample", function(object, value) standardGeneric("set_info_sample"))

setMethod("set_info_sample", "mgnet", function(object, value) {
  object@info_sample <- value
  validObject(object)
  return(object)
})

setMethod("set_info_sample", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_info_sample(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})


# SET LINEAGE
#------------------------------------------------------------------------------#
#' Set `lineage` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the lineage data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new lineage data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated lineage data.
#' @export
#' @name set_lineage
#' @aliases set_lineage,mgnet-method set_lineage,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_lineage", function(object, value) standardGeneric("set_lineage"))

setMethod("set_lineage", "mgnet", function(object, value) {
  object@lineage <- value
  validObject(object)
  return(object)
})

setMethod("set_lineage", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_lineage(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})


# SET INFO_TAXA
#------------------------------------------------------------------------------#
#' Set `info_taxa` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the info_taxa data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new info_taxa data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated info_taxa data.
#' @export
#' @name set_info_taxa
#' @aliases set_info_taxa,mgnet-method set_info_taxa,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_info_taxa", function(object, value) standardGeneric("set_info_taxa"))

setMethod("set_info_taxa", "mgnet", function(object, value) {
  object@info_taxa <- value
  validObject(object)
  return(object)
})

setMethod("set_info_taxa", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_info_taxa(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})



# SET NETWORK
#------------------------------------------------------------------------------#
#' Set `network` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the network data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new network data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated network data.
#' @export
#' @name set_network
#' @aliases set_network,mgnet-method set_network,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_network", function(object, value) standardGeneric("set_network"))

setMethod("set_network", "mgnet", function(object, value) {
  object@network <- value
  validObject(object)
  return(object)
})

setMethod("set_network", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_network(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})


# SET COMMUNITY
#------------------------------------------------------------------------------#
#' Set `community` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the community data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new community data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated community data.
#' @export
#' @name set_community
#' @aliases set_community,mgnet-method set_community,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_community", function(object, value) standardGeneric("set_community"))

setMethod("set_community", "mgnet", function(object, value) {
  object@community <- value
  validObject(object)
  return(object)
})

setMethod("set_community", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_community(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})