setOldClass("igraph")
setOldClass("communities")
################################################################################
################################################################################
# CLASS MGNET
################################################################################
################################################################################
#' S4 class to manage metagenomic networks
#'
#'
#'
#' @slot data table of data
#' @slot meta sample features
#' @slot taxa taxonomic information
#' @slot netw undirected, weighted and signed igraph network
#' @slot comm graph communities of netw
#'
#' @import methods
#' @importFrom igraph make_empty_graph cluster_fast_greedy V vcount
#' @name mgnet-class
#' @rdname mgnet-class
#' @exportClass mgnet
mgnet <- setClass(
  Class="mgnet",
  
  contains="mg",
  
  slot=c(
    netw="igraph",
    comm="communities"),
  
  prototype=prototype(data=matrix(nrow=0,ncol=0),
                      meta=data.frame(),
                      taxa=matrix(nrow=0,ncol=0),
                      netw=make_empty_graph(n=0, directed=FALSE),
                      comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))),
  
  validity=function(object){
    
    # CHECK NETW
    #-------------------------------------#
    if(length(object@netw)!=0){
      if(is_directed(object@netw)) return("\n netw slot must be an undirected igraph object")
      if(!is_weighted(object@netw)) return("\n netw slot must be an weighted igraph object")
    }
    
    # CHECK NETW AND DATA
    #-------------------------------------#
    if(length(object@netw)!=0 & length(object@data)!=0){
      if(vcount(object@netw)!=ncol(object@data)) return("\n netw vertices number must be equal to columns number of data")
      if(any(colnames(object@data)!=V(object@netw)$name)) return("\n netw vertex name and colnames of data must be equal")
    }
    
    # CHECK NETW AND TAXA
    #-------------------------------------#
    if(length(object@netw)!=0 & length(object@taxa)!=0){
      if(vcount(object@netw)!=nrow(object@taxa)) return("\n netw vertices number must be equal to rows number of taxa")
      if(any(rownames(object@taxa)!=V(object@netw)$name)) return("\n netw vertex name and rownames of taxa must be equal")
    }
    
    # CHECK COMM
    #-------------------------------------#
    if(length(object@comm)!=0 & length(object@netw)==0) return("\n you can't have the comm slot without its associated netw")
    if(length(object@comm)!=0){
      if(object@comm$vcount!=vcount(object@netw)) stop("comm and netw must have the same number of vertices")
    }
    
    TRUE
  }
)




################################################################################
################################################################################
# CONSTRUCTOR MGNET
################################################################################
################################################################################
#' User constructor for mgnet s4 class.
#' 
#' @rdname mgnet-constructor
#' @aliases mgnet,constructor
#' 
#' @description User constructor to create an object belonging to formal s4 
#' class mgnet avoiding the new function.
#' 
#' @param data numeric matrix with all elements >=0.  
#' @param meta data.frame with experimental variables.
#' @param taxa character matrix with taxonomic classification.  
#' @param netw un-directed weighted igraph network.
#' @param comm communities object associated to netw.
#' @param mg object belong to mg class.
#' @export
mgnet <- function(data=matrix(nrow=0,ncol=0),
                  meta=data.frame(),
                  taxa=matrix(nrow=0,ncol=0),
                  netw=make_empty_graph(n=0, directed=FALSE),
                  comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE)),
                  mg=new("mg")){
  
  if((length(data)!=0 | length(meta)!=0 | length(taxa)!=0) & !empty(mg)){
    stop("if using an object of class mg then data, meta and taxa must be left empty")
  }
  
  if(empty(mg)){
    return(new("mgnet",data=data,meta=meta,taxa=taxa,netw=netw,comm=comm))
  } else {
    return(new("mgnet",data=mg@data,meta=mg@meta,taxa=mg@taxa,netw=netw,comm=comm))
  }
  
}
################################################################################
################################################################################
# END CONSTRUCTOR MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# GETTERS MGNET
################################################################################
################################################################################
# NETW
#####################################
#' Retrieves netw.
#' 
#' @description 
#' Return the igraph un-directed weighted network.
#'
#' @usage netw(object)
#'
#' @param object (Required) \code{\link{mgnet-class}}.
#'
#' @rdname netw
#' @docType methods
#' @export
setGeneric("netw", function(object) standardGeneric("netw"))
#' @rdname netw
#' @aliases netw,mgnet
setMethod("netw", "mgnet", function(object){return(object@netw)})
#####################################
# COMM
#####################################
#' Retrieves network communities.
#'
#' @description
#' Return the communities associated to netw slot.
#'
#' @usage comm(object)
#'
#' @param object (Required) \code{\link{mgnet-class}}.
#'
#' @rdname comm
#' @docType methods
#' @export
setGeneric("comm", function(object) standardGeneric("comm"))
#' @rdname comm
#' @aliases comm,mgnet
setMethod("comm", "mgnet", function(object){return(object@comm)})
################################################################################
################################################################################
# END GETTERS MGNET
################################################################################
################################################################################



################################################################################
################################################################################
# SETTERS MGNET
################################################################################
################################################################################
# DATA<-
#####################################
#' @rdname assign-data
#' @aliases data<-,mgnet,matrix
setMethod("data<-", c("mgnet", "matrix"), function(object, value){
  new("mgnet",data=value, meta=object@meta, taxa=object@taxa, netw=object@netw,
      comm=object@comm)
})
#####################################
# META<-
#####################################
#' @rdname assign-meta
#' @aliases meta<-,mgnet,data.frame
setMethod("meta<-", c("mgnet", "data.frame"), function(object, value){
  new("mgnet",data=object@data, meta=value, taxa=object@taxa, netw=object@netw,
      comm=object@comm)
})
#####################################
# TAXA<-
#####################################
#' @rdname assign-taxa
#' @aliases taxa<-,mgnet,matrix
setMethod("taxa<-", c("mgnet", "matrix"), function(object, value){
  new("mgnet",data=object@data, meta=object@meta, taxa=value, netw=object@netw,
      comm=object@comm)
})
#####################################
# NETW<-
#####################################
#' Assign a new netw to \code{object}
#'
#' @usage netw(object) <- value
#'
#' @param object (Required) \code{\link{mgnet-class}}.
#' @param value (Required) \code{\link{igraph}}
#'
#' @export
#' @docType methods
#' @rdname assign-netw
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))
#' @rdname assign-netw
#' @aliases netw<-,mgnet,igraph
setMethod("netw<-", c("mg", "igraph"), function(object, value){
  new("mgnet",data=object@data, meta=object@meta, taxa=object@taxa, netw=value,
      comm=object@comm)
})
#####################################
# COMM<-
#####################################
#' Assign a new comm to \code{object}
#'
#' @usage comm(object) <- value
#'
#' @param object (Required) \code{\link{mgnet-class}}.
#' @param value (Required) \code{\link{communities}}
#'
#' @export
#' @docType methods
#' @rdname assign-comm
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))
#' @rdname assign-comm
#' @aliases comm<-,mgnet,communities
setMethod("comm<-", c("mg", "communities"), function(object, value){
  new("mgnet",data=object@data, meta=object@meta, taxa=object@taxa,
      netw=object@netw, comm=value)
})
################################################################################
################################################################################
# END SETTERS MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# EXTRACTOR MGNET
################################################################################
################################################################################
#' Method extensions to extraction operator for mgnet object.
#'
#' @param x See \code{\link[base]{Extract}}, \code{\link{mg}} object.
#' @param i See \code{\link[base]{Extract}}, samples indices.
#' @param j See \code{\link[base]{Extract}}, taxa indices
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
setMethod(f="[",
          signature="mg",
          definition=function(x,i,j){
            warning("operator '[' is not well defined in mgnet class.
                    It return a simply mg without netw and comm info")
            return(new("mg",
                       data=x@data[i,j,drop=FALSE],
                       meta=x@meta[i, ,drop=FALSE],
                       taxa=x@taxa[j, ,drop=FALSE]))
          })
################################################################################
################################################################################
# END EXTRACTOR MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# SHOW METHOD MGNET
################################################################################
################################################################################
#'@importFrom igraph ecount edge_density membership sizes
setMethod("show","mg",
          function(object){
            cat("******* Class mgnet , method Show ******* \n")
            cat(paste("Sample Number:",max(nrow(object@data),nrow(object@meta)),"\n"))
            cat(paste("Taxa Number:",max(ncol(object@data),nrow(object@taxa)),"\n"))
            cat(paste("Sample Meta Data:",paste(colnames(object@meta),collapse="," )),"\n")
            cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")),"\n")
            cat(paste("Edge Number (Density): ",ecount(object@netw)," (",
                      round(edge_density(object@netw)*100,2),"%)","\n",sep=""))
            
            if(length(object@comm)!=0){
              cat(paste("Signed Communities Number:",max(membership(object@comm))))
              
              if("0" %in% names(sizes(object@comm))){
                cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=",")))
                cat(paste("Isolated Nodes:", sizes(object@comm)[[1]]))
              } else {
                cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=",")))
                cat("There aren't isolated nodes")
              }
            } else {
              cat("There aren't communities information \n")
            }
            
            cat("********** End Show (mgnet) ********** \n")
          })
################################################################################
################################################################################
# END SHOW METHOD MGNET
################################################################################
################################################################################



