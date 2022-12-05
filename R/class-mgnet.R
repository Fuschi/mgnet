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