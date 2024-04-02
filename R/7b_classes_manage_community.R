# REMOVE SMALLER COMMUNITIES
#------------------------------------------------------------------------------#
#' Remove or Isolate Smaller Communities
#'
#' This function processes an object (either `igraph`, `mgnet`, or `mgnetList`) to remove or isolate communities based on their size,
#' with additional control over the output format regarding how isolated nodes are labeled.
#'
#' @description
#' The function `remove_smaller_comm` removes or isolates vertices from communities that are below a specified size threshold.
#' It allows further customization of the output by optionally labeling isolated nodes in a specified format.
#'
#' @param object An object of class `igraph`, `mgnet`, or `mgnetList` representing the network.
#' @param size A positive integer indicating the minimum number of vertices required for a community to be retained.
#' @param trim A logical value indicating whether vertices from smaller communities should be removed (`TRUE`)
#'        or simply isolated within the network (`FALSE`).
#' @param membership_format A character string specifying the format of the output regarding isolated nodes.
#'        `"mgnet"` will label isolated nodes as `"0"`, whereas `"igraph"` retains the default behavior
#'        where isolated nodes are not specifically labeled. This parameter is particularly relevant when `trim` is `FALSE`.
#' @return
#' The function returns an object of the same class as the input (`object`) but with smaller communities either removed or isolated
#' based on the specified `size` and `trim` parameters. The output format of isolated nodes can be adjusted with the `membership_format` parameter.
#'
#' @export
#' @aliases remove_smaller_communities,mgnet-method
#' @importFrom igraph delete_vertices subgraph.edges
setGeneric("remove_smaller_communities", function(object, size, trim = TRUE, membership_format = "mgnet") standardGeneric("remove_smaller_communities"))

setMethod("remove_smaller_communities", "mgnet", function(object, size,
                                                          trim = TRUE, membership_format = "mgnet"){

  #Checks Arguments
  if( length(object@network)==0 || length(object@community)==0 ) stop("network and community slots must be present")
  if( !is.numeric(size) || round(size) != size || size < 0 ) stop("size must be integer greater than 0")
  if( !is.logical(trim) ) stop("keep must be logical")
  if( length(object)==0 ) stop("missing communities info")
  membership_format <- match.arg(membership_format, c("mgnet","igraph"))
  #End Checks

  graph <- object@network
  comm <- as_mgnet_communities(object@community)

  keep.comm.names <- as.numeric(names(sizes(comm)[sizes(comm)>=size]))
  keep.comm.vids <-  which(comm$membership %in% keep.comm.names)

  if(!trim){
    
    community(object)$membership[setdiff(1:comm$vcount, keep.comm.vids)] <- 0
    community(object)$modularity <- NA
    if( membership_format == "igraph") community(object) <- as_igraph_communities(community(object))
    validObject(object)
    return(object)
    
  } else {
    
    object <- object[ ,keep.comm.vids]
    validObject(object)
    return(object)
  }
  
})

# TO UPDATE
#' # ISOLATE 
#' #------------------------------------------------------------------------------#
#' #' Isolate Vertices
#' #' 
#' #' @description Remove vertices from their communities.
#' #' 
#' #' @usage isolate(object, ..., condition="AND")
#' #' 
#' #' @param object mgnet.
#' #' @param ... arbitrary set of vectors of functions. If vectors they must contain
#' #' numeric or name or logical indices to be preserved. If functions they must return
#' #' indices as just explained.
#' #' @param condition (Optional default "AND") character indicates the logical operation between indices.
#' #' 
#' #' @importFrom igraph make_clusters
#' #' @rdname isolate-methods
#' #' @docType methods
#' #' @export
#' setGeneric("isolate", function(object,...,condition="AND") standardGeneric("isolate"))
#' #' @rdname isolate-methods
#' setMethod("isolate", "mgnet", function(object,...,condition="AND"){
#'   
#'   condition <- match.arg(condition, c("AND","OR"))
#'   if(length(object@netw)==0) stop("netw slot cannot be empty")
#'   
#'   
#'   IDX <- list(...)
#'   for(i in 1:length(IDX)){
#'     
#'     idx <- IDX[[i]]
#'     if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
#'     
#'     if(is.function(idx)){
#'       idx <- idx(object)
#'       if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
#'     }  
#'     
#'     
#'     if(is.numeric(idx)){
#'       #numeric
#'       if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>ntaxa(object)){
#'         stop("numeric indices must be integers in range 1 to the maximum nuber of taxa")
#'       }
#'       #transform to logical
#'       IDX[[i]] <- (1:nsample(object))%in%idx
#'       
#'     }else if(is.character(idx)){
#'       #character
#'       if( !any(idx%in%taxa_name(object)) ){
#'         stop("string indices must be a subset of taxa_name of object")
#'       }
#'       #transform to logical
#'       IDX[[i]] <- taxa_name(object)%in%idx
#'       
#'     }else if(is.logical(idx)){
#'       #logical
#'       if(length(idx)!=ntaxa(object)){
#'         stop("logical indices must have the length equal to taxa number")
#'       }
#'       IDX[[i]] <- idx
#'       
#'     } else {stop("the indices must character, numeric or logical vectors")}
#'   }
#'   
#'   
#'   IDX <- do.call(rbind,IDX)
#'   # Condition  
#'   if(condition=="AND"){
#'     IDX <- as.logical(apply(IDX,2,prod))
#'   } else {
#'     IDX <- apply(IDX,2,sum)>0
#'   }
#'   
#'   membership <- comm(object)$membership
#'   membership[IDX] <- 0
#'   
#'   comm.new <- make_clusters(netw(object),
#'                             membership=membership,
#'                             algorithm="modified",
#'                             modularity=1)
#'   comm.new$modularity <- NA
#'   
#'   return(mgnet(data=object@data,
#'                meta_sample=object@meta_sample,
#'                taxa=object@taxa,
#'                meta_taxa=object@meta_taxa,
#'                log_data=object@log_data,
#'                netw=object@netw,
#'                comm=comm.new))
#'   
#' })
#' #' @rdname isolate-methods
#' setMethod("isolate","list",
#'           function(object,...,condition="AND"){
#'             is_mgnet_list(object)
#'             lapply(object, selectMethod(f="isolate",
#'                                         signature="mgnet"),
#'                    ...=..., condition=condition)}
#' )
#' ################################################################################
#' ################################################################################
#' # END ISOLATE
#' ################################################################################
#' ################################################################################