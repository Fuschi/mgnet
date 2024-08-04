#' Mutate Sample Information based on Abundance Criteria
#'
#' This function allows for mutating sample information based on specified criteria applied to
#' the abundance, relative abundance, or normalized abundance data in an `mgnet` object.
#' Custom functions define the criteria for mutation, enabling flexible modification of sample
#' information across multiple datasets.
#'
#' @param object An `mgnetList` object.
#' @param abundance A named list of functions to be applied to the abundance data matrix.
#'        Each function should take a single row of abundance data as input and return the value
#'        to be stored in the `info_sample` slot. Default is `NULL`.
#' @param rel_abundance A named list of functions to be applied to the relative abundance data matrix.
#'        Each function should take a single row of relative abundance data as input and return
#'        the value to be stored in the `info_sample` slot. Default is `NULL`.
#' @param norm_abundance A named list of functions to be applied to the log-transformed abundance data
#'        matrix. Each function should take a single row of log-transformed abundance data as input
#'        and return the value to be stored in the `info_sample` slot. Default is `NULL`.
#' @param mgnet A named list of functions to be applied direct on the mgnet object.
#'
#' @return A modified `mgnet` object with the `info_sample` slot updated based on the applied mutations.
#'
#' @export
#' @seealso \link{filter_criteria_sample} for sample-based filtering using criteria.
#' @aliases mutate_criteria_sample,mgnet-method mutate_criteria_sample,mgnetList-method
setGeneric("mutate_criteria_sample",
           function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                    mgnet = NULL) {
             standardGeneric("mutate_criteria_sample")
           })

setMethod("mutate_criteria_sample", "mgnet",
          function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                   mgnet = NULL) {
            
            # Checks
            if(is.null(mgnet) & is.null(abundance) & is.null(rel_abundance) & is.null(norm_abundance)) return(object)
            
            # Check if all mutations are named list of functions or NULL
            check_mutations <- function(mutations) {
              if (is.null(mutations)) return(TRUE)
              if (!is.list(mutations)) return(FALSE)
              if (!all(sapply(names(mutations), nzchar))) return(FALSE)
              if (!all(sapply(mutations, is.function))) return(FALSE)
              TRUE
            }
            if (!all(sapply(list(mgnet, abundance, rel_abundance, norm_abundance), check_mutations))) {
              stop("All criteria must be named list of functions or NULL.")
            }
            
            apply_mutations <- function(data, mutations){
              
              if(length(data)==0) stop(paste(data, "slot missing."))
              if (is.null(mutations)) return(data.frame())
              
              result <- sapply(mutations, \(f){
                apply(data,1,f)
              }, USE.NAMES = TRUE, simplify = TRUE)
              
              if((!is.matrix(result) | nrow(result) != nsample(object) | ncol(result) != length(mutations))){
                stop(paste("All functions associated to", data, "must return a single vector with length equal to the sample number of the object"))
              }
              return(result)
            }
            
            apply_mutations_mgnet <- function(data, mutations){
              
              if(length(data)==0) stop(paste(data, "slot missing."))
              if (is.null(mutations)) return(data.frame())
              
              result <- sapply(mutations, \(f){
                f(data)
              }, USE.NAMES = TRUE, simplify = TRUE)
              
              if((!is.matrix(result) | nrow(result) != nsample(object) | ncol(result) != length(mutations))){
                stop(paste("All functions associated to", data, "must return a single vector with length equal to the sample number of the object"))
              }
              return(result)
            }
            
            
            
            # Apply mutations functions to each data type
            mgnet_new <- if(length(object)!=0){
              apply_mutations_mgnet(object, mgnet)
            } else {
              data.frame()
            }
            abundance_new <- if(length(object@abundance)!=0){
              apply_mutations(object@abundance, abundance)
            } else {
              data.frame()
            }
            rel_abundance_new <- if(length(object@rel_abundance)!=0){
              apply_mutations(object@rel_abundance, rel_abundance)
            } else {
              data.frame()
            }
            norm_abundance_new <- if(length(object@norm_abundance)!=0){
              apply_mutations(object@norm_abundance, norm_abundance)
            } else {
              data.frame()
            }
            
            missing_cbind <- function(df1, df2){
              
              if(length(df1)==0 && length(df2)==0) return(data.frame())
              if(length(df1)==0 && length(df2)!=0) return(df2)
              if(length(df1)!=0 && length(df2)==0) return(df1)
              if(length(df1)!=0 && length(df2)!=0) return(cbind(df1, df2))
              
            }
            
            result <- missing_cbind(abundance_new, rel_abundance_new)
            result <- missing_cbind(result, norm_abundance_new)
            result <- missing_cbind(result, mgnet_new)
            result <- missing_cbind(info_sample(object), result)
            
            info_sample(object) <- as.data.frame(result)
            
            return(object)
          })

setMethod("mutate_criteria_sample", "mgnetList",
          function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                   mgnet = NULL) {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              mutate_criteria_sample(mgnet_obj, abundance, rel_abundance, norm_abundance,
                                     mgnet)
            })
            return(object)
          })



#' Mutate Taxa Information based on Abundance Criteria
#'
#' This function allows for mutating taxa information based on specified criteria applied to
#' the abundances matrices in an `mgnet` object or `mgnetList` object.
#' Custom functions define the criteria for mutation, enabling flexible modification of taxa
#' information across multiple datasets.
#'
#' @param object An `mgnetList` object.
#' @param abundance A named list of functions to be applied to the abundance data matrix.
#'        Each function should take a single row of abundance data as input and return the value
#'        to be stored in the `info_taxa` slot. Default is `NULL`.
#' @param rel_abundance A named list of functions to be applied to the relative abundance data matrix.
#'        Each function should take a single row of relative abundance data as input and return
#'        the value to be stored in the `info_taxa` slot. Default is `NULL`.
#' @param norm_abundance A named list of functions to be applied to the log-transformed abundance data
#'        matrix. Each function should take a single row of log-transformed abundance data as input
#'        and return the value to be stored in the `info_taxa` slot. Default is `NULL`.
#' @param network A named list of functions to be applied to the network data matrix.
#'        Each function should take a single column of network data as input and return
#'        the value to be stored in the `info_taxa` slot. Default is `NULL`.
#' @param community A named list of functions to be applied to the community data matrix.
#'        Each function should take a single column of community data as input and return
#'        the value to be stored in the `info_taxa` slot. Default is `NULL`.
#' @param mgnet A named list of functions to be applied direct on the mgnet object.
#'
#' @return A modified `mgnet` object with the `info_taxa` slot updated based on the applied mutations.
#'
#' @export
#' @seealso \link{filter_criteria_taxa} for taxa-based filtering using criteria.
#' @aliases mutate_criteria_taxa,mgnet-method mutate_criteria_taxa,mgnetList-method
setGeneric("mutate_criteria_taxa",
           function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                    network = NULL, community = NULL, mgnet = NULL) {
             standardGeneric("mutate_criteria_taxa")
           })

setMethod("mutate_criteria_taxa", "mgnet",
          function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL, 
                   network = NULL, community = NULL, mgnet = NULL) {
            
            # Checks
            if(is.null(abundance) & is.null(rel_abundance) & is.null(norm_abundance) & is.null(community) & is.null(network) & is.null(mgnet)) return(object)
            
            # Check if all mutations are named list of functions or NULL
            check_mutations <- function(mutations) {
              if (is.null(mutations)) return(TRUE)
              if (!is.list(mutations)) return(FALSE)
              if (!all(sapply(names(mutations), nzchar))) return(FALSE)
              if (!all(sapply(mutations, is.function))) return(FALSE)
              TRUE
            }
            if (!all(sapply(list(mgnet, abundance, rel_abundance, norm_abundance, community, network), check_mutations))) {
              stop("All criteria must be named list of functions or NULL.")
            }
            
            apply_mutations <- function(data, mutations){
              
              if(length(data)==0) stop(paste(data, "slot missing."))
              if (is.null(mutations)) return(data.frame())
              
              result <- sapply(mutations, \(f){
                apply(data,2,f)
              }, USE.NAMES = TRUE, simplify = TRUE)
              
              if((!is.matrix(result) | nrow(result) != ntaxa(object) | ncol(result) != length(mutations))){
                stop(paste("All functions associated to", data, "must return a single vector with length equal to the taxa number of the object"))
              }
              return(result)
            }
            
            apply_mutations_mgnet_netw_comm <- function(data, mutations){
              
              if(length(data)==0) stop(paste(data, "slot missing."))
              if (is.null(mutations)) return(data.frame())
              
              result <- sapply(mutations, \(f){
                f(data)
              }, USE.NAMES = TRUE, simplify = TRUE)
              
              if((!is.matrix(result) | nrow(result) != ntaxa(object) | ncol(result) != length(mutations))){
                stop(paste("All functions must return a single vector with length equal to the sample number of the object"))
              }
              return(result)
            }
            
            # Apply mutations functions to each data type
            abundance_new <- if(length(object@abundance)!=0){
              apply_mutations(object@abundance, abundance)
            } else {
              data.frame()
            }
            rel_abundance_new <- if(length(object@rel_abundance)!=0){
              apply_mutations(object@rel_abundance, rel_abundance)
            } else {
              data.frame()
            }
            norm_abundance_new <- if(length(object@norm_abundance)!=0){
              apply_mutations(object@norm_abundance, norm_abundance)
            } else {
              data.frame()
            }
            network_new <- if(length(object@network)!=0){
              apply_mutations_mgnet_netw_comm(object@network, network)
            } else {
              data.frame()
            }
            community_new <- if(length(object@community)!=0){
              apply_mutations_mgnet_netw_comm(object@community, community)
            } else {
              data.frame()
            }
            mgnet_new <- if(length(object)!=0){
              apply_mutations_mgnet_netw_comm(object, mgnet)
            } else {
              data.frame()
            }
            
            missing_cbind <- function(df1, df2){
              
              if(length(df1)==0 && length(df2)==0) return(data.frame())
              if(length(df1)==0 && length(df2)!=0) return(df2)
              if(length(df1)!=0 && length(df2)==0) return(df1)
              if(length(df1)!=0 && length(df2)!=0) return(cbind(df1, df2))
              
            }
            
            result <- missing_cbind(abundance_new, rel_abundance_new)
            result <- missing_cbind(result, norm_abundance_new)
            result <- missing_cbind(result, network_new)
            result <- missing_cbind(result, community_new)
            result <- missing_cbind(result, mgnet_new)
            result <- missing_cbind(info_taxa(object), result)
            
            
            info_taxa(object) <- as.data.frame(result)
            
            return(object)
          })

setMethod("mutate_criteria_taxa", "mgnetList",
          function(object, abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                   network = NULL, community = NULL, mgnet = NULL) {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              mutate_criteria_taxa(mgnet_obj, abundance, rel_abundance, norm_abundance, network, community, mgnet)
            })
            return(object)
          })