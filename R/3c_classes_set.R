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


# SET REL_ABUNDANCE
#------------------------------------------------------------------------------#
#' Set `rel_abundance` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the rel_abundance data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new rel_abundance data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated rel_abundance data.
#' @export
#' @name set_rel_abundance
#' @aliases set_rel_abundance,mgnet-method set_rel_abundance,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_rel_abundance", function(object, value) standardGeneric("set_rel_abundance"))

setMethod("set_rel_abundance", "mgnet", function(object, value) {
  object@rel_abundance <- value
  validObject(object)
  return(object)
})

setMethod("set_rel_abundance", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_rel_abundance(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})

# SET norm_abundance
#------------------------------------------------------------------------------#
#' Set `norm_abundance` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the rel_abundance data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new norm_abundance data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated rel_abundance data.
#' @export
#' @name set_norm_abundance
#' @aliases set_norm_abundance,mgnet-method set_norm_abundance,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_norm_abundance", function(object, value) standardGeneric("set_norm_abundance"))

setMethod("set_norm_abundance", "mgnet", function(object, value) {
  object@norm_abundance <- value
  validObject(object)
  return(object)
})

setMethod("set_norm_abundance", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in seq_along(object@mgnets)) {
    object@mgnets[[i]] <- set_norm_abundance(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})

# SET CLR_abundance
#------------------------------------------------------------------------------#
#' Store clr-transformed data in `norm_abundance` Slot in `mgnet` Objects
#'
#' Applies centered log-ratio (CLR) transformation or inter-quantile log-ratio (ICLR) transformation to the abundance data
#' in `mgnet` objects, after handling zeros with specified strategies. This transformation is useful for compositional data analysis,
#' making the data suitable for statistical modeling and comparison.
#'
#' @param object An `mgnet` object.
#' @param clr_variant The method to calculate norm_abundance. Options are 'clr' for centered log-ratio and 'iclr' for
#'        inter-quantile log-ratio. Default is 'clr'.
#' @param zero_strategy The zero replacement strategy before log-ratio transformation.
#'        'const' replaces zeros with a small constant, and 'unif' replaces zeros with a uniform random value
#'        between 6.5% and 65% of the detection limit. Default is 'unif'.
#' @return The `mgnet` object with updated `norm_abundance` data.
#' @export
#' @seealso \code{\link{clr}}, \code{\link{iclr}}, \code{\link{zero_dealing}}
#'         See also:
#'         \code{\link[=zero_dealing]{zero_dealing}} for zero replacement strategies,
#'         \code{\link[=clr]{clr}} for details on centered log-ratio transformation,
#'         \code{\link[=iclr]{iclr}} for details on inter-quantile log-ratio transformation.
#'         
#' @name set_CLR_abundance
#' @aliases set_CLR_abundance,mgnet-method set_CLR_abundance,mgnetList-method
setGeneric("set_CLR_abundance", function(object, clr_variant = "clr", zero_strategy = "unif") standardGeneric("set_CLR_abundance"))

setMethod("set_CLR_abundance", "mgnet", function(object, clr_variant = "clr", zero_strategy = "unif") {
  
  clr_variant <- match.arg(clr_variant, c("clr", "iclr"))
  zero_strategy <- match.arg(zero_strategy, c("unif", "const"))
  if(length(object@abundance) == 0) stop("Abundance matrix missing; cannot calculate log-ratio abundance.")
  
  abundance_nozero <- zero_dealing(object@abundance, method = zero_strategy)
  
  if(clr_variant == "clr") {
    object@norm_abundance <- clr(abundance_nozero)
  } else {
    object@norm_abundance <- iclr(abundance_nozero)
  }
  
  validObject(object)
  return(object)
})

setMethod("set_CLR_abundance", "mgnetList", function(object, clr_variant = "clr", zero_strategy = "unif") {
  object <- lapply(object@mgnets, function(mgnet) {
    set_CLR_abundance(mgnet, clr_variant, zero_strategy)
  })
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