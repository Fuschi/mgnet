# SET ABUNDANCE
#------------------------------------------------------------------------------#
#' Set `abun` Slot in `mgnet` Objects
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
#' @name set_abun
#' @aliases set_abun,mgnet-method set_abun,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_abun", function(object, value) standardGeneric("set_abun"))

setMethod("set_abun", "mgnet", function(object, value) {
  abun(object) <- value
  return(object)
})

setMethod("set_abun", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  abun(object) <- value
  return(object)
})


# SET rela
#------------------------------------------------------------------------------#
#' Set `rela` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the rela data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new rela data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated rela data.
#' @export
#' @name set_rela
#' @aliases set_rela,mgnet-method set_rela,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_rela", function(object, value) standardGeneric("set_rela"))

setMethod("set_rela", "mgnet", function(object, value) {
  rela(object) <- value
  return(object)
})

setMethod("set_rela", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  rela(object) <- value
  return(object)
})

# SET norm
#------------------------------------------------------------------------------#
#' Set `norm` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the rela data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new norm data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated rela data.
#' @export
#' @name set_norm
#' @aliases set_norm,mgnet-method set_norm,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_norm", function(object, value) standardGeneric("set_norm"))

setMethod("set_norm", "mgnet", function(object, value) {
  norm(object) <- value
  return(object)
})

setMethod("set_norm", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  norm(object) <- value
  return(object)
})

# SET NORM CLR
#------------------------------------------------------------------------------#
#' Store clr-transformed data in `norm` Slot in `mgnet` Objects
#'
#' Applies centered log-ratio (CLR) transformation or inter-quantile log-ratio (ICLR) transformation to the abundance data
#' in `mgnet` objects, after handling zeros with specified strategies. This transformation is useful for compositional data analysis,
#' making the data suitable for statistical modeling and comparison.
#'
#' @param object An `mgnet` object.
#' @param clr_variant The method to calculate norm. Options are 'clr' for centered log-ratio and 'iclr' for
#'        inter-quantile log-ratio. Default is 'clr'.
#' @param zero_strategy The zero replacement strategy before log-ratio transformation.
#'        'const' replaces zeros with a small constant, and 'unif' replaces zeros with a uniform random value
#'        between 6.5% and 65% of the detection limit. Default is 'const'.
#' @return The `mgnet` object with updated `norm` data.
#' @export
#' @seealso \code{\link{clr}}, \code{\link{iclr}}, \code{\link{zero_dealing}}
#'         See also:
#'         \code{\link[=zero_dealing]{zero_dealing}} for zero replacement strategies,
#'         \code{\link[=clr]{clr}} for details on centered log-ratio transformation,
#'         \code{\link[=iclr]{iclr}} for details on inter-quantile log-ratio transformation.
#'         
#' @name set_norm_CLR
#' @aliases set_norm_CLR,mgnet-method set_norm_CLR,mgnetList-method
setGeneric("set_norm_CLR", function(object, clr_variant = "clr", zero_strategy = "const") standardGeneric("set_norm_CLR"))

setMethod("set_norm_CLR", "mgnet", function(object, clr_variant = "clr", zero_strategy = "const") {
  
  clr_variant <- match.arg(clr_variant, c("clr", "iclr"))
  zero_strategy <- match.arg(zero_strategy, c("const", "unif"))
  if(length(object@abun) == 0 & length(object@rela)) stop("abundance and relative data missing; cannot calculate log-ratio abundance.")
  
  if(length(object@abun) != 0){
    abundance <- abun(object, .fmt = "mat")
  } else {
    abundance <- rela(object, .fmt = "mat")
  }
  
  norm(object) <- clr_zero_handle(X = abundance, clr_variant = clr_variant, zero_strategy = zero_strategy)
  return(object)
})

setMethod("set_norm_CLR", "mgnetList", function(object, clr_variant = "clr", zero_strategy = "const") {
  norm(object) <- sapply(object, \(x){
    
    clr_variant <- match.arg(clr_variant, c("clr", "iclr"))
    zero_strategy <- match.arg(zero_strategy, c("const", "unif"))
    if(length(x@abun) == 0 & length(x@rela)) stop("abundance and relative data missing; cannot calculate log-ratio abundance.")
    
    if(length(x) != 0){
      abundance <- abun(x, .fmt = "mat")
    } else {
      abundance <- rela(x, .fmt = "mat")
    }
    
    return(clr_zero_handle(X = abundance, clr_variant = clr_variant, zero_strategy = zero_strategy))
  }, simplify = FALSE, USE.NAMES = TRUE)
  return(object)
})

# SET META
#------------------------------------------------------------------------------#
#' Set `meta` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the sample metadata for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new sample metadata data to be set, a data.frame for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated sample data.
#' @export
#' @name set_meta
#' @aliases set_meta,mgnet-method set_meta,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_meta", function(object, value) standardGeneric("set_meta"))

setMethod("set_meta", "mgnet", function(object, value) {
  object@sample <- value
  validObject(object)
  return(object)
})

setMethod("set_meta", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  
  for (i in names(object)) {
    object@mgnets[[i]] <- set_meta(object@mgnets[[i]], value[[i]])
    validObject(object@mgnets[[i]])
  }
  
  validObject(object)
  return(object)
})


# SET taxa
#------------------------------------------------------------------------------#
#' Set `taxa` Slot in `mgnet` Objects
#'
#' @description
#' This function sets the taxa data for an `mgnet` object or each `mgnet` object 
#' within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param value The new taxa data to be set, a numeric matrix for `mgnet` objects 
#' or a list of numeric matrices for `mgnetList` objects.
#' @return The modified `mgnet` or `mgnetList` object with the updated taxa data.
#' @export
#' @name set_taxa
#' @aliases set_taxa,mgnet-method set_taxa,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_taxa", function(object, value) standardGeneric("set_taxa"))

setMethod("set_taxa", "mgnet", function(object, value) {
  taxa(object) <- value
  return(object)
})

setMethod("set_taxa", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  taxa(object) <- value
  return(object)
})


# SET netw
#------------------------------------------------------------------------------#
#' Set `netw` Slot in `mgnet` Objects
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
#' @name set_netw
#' @aliases set_netw,mgnet-method set_netw,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_netw", function(object, value) standardGeneric("set_netw"))

setMethod("set_netw", "mgnet", function(object, value) {
  netw(object) <- value
  return(object)
})

setMethod("set_netw", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  netw(object) <- value
  return(object)
})


# SET COMMUNITY
#------------------------------------------------------------------------------#
#' Set `comm` Slot in `mgnet` Objects
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
#' @name set_comm
#' @aliases set_comm,mgnet-method set_comm,mgnetList-method
#' @importFrom methods validObject
setGeneric("set_comm", function(object, value) standardGeneric("set_comm"))

setMethod("set_comm", "mgnet", function(object, value) {
  comm(object) <- value
  return(object)
})

setMethod("set_comm", "mgnetList", function(object, value) {
  are_list_assign(object, value)
  comm(object) <- value
  return(object)
})
