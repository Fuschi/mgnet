setOldClass("igraph")
setOldClass("communities")
################################################################################
################################################################################
# CLASS MGNET
################################################################################
################################################################################
#' S4 class to manage metagenomics networks
#'
#' @slot data ngs data
#' @slot meta_sample samples experimental features.
#' @slot taxa taxonomic table 
#' @slot meta_taxa additional info on taxa.
#' @slot netw undirected, weighted igraph network.
#' @slot comm graph communities of netw.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy V vcount
#' @import methods
#' @name mgnet-class
#' @rdname mgnet-class
#' @exportClass mgnet
setClass(
  Class="mgnet",
  
  slot=c(data="matrix",
         meta_sample="data.frame",
         taxa="matrix",
         meta_taxa="data.frame",
         netw="igraph",
         comm="communities"
         ),
  
  prototype=prototype(data=matrix(nrow=0,ncol=0),
                      meta_sample=data.frame(),
                      taxa=matrix(nrow=0,ncol=0),
                      meta_taxa=data.frame(),
                      netw=make_empty_graph(n=0, directed=FALSE),
                      comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))),
  
  validity=function(object){
    
    #CHECK DATA
    #-------------------------------------#
    if( length(object@data)!=0 ){
      if(!is.numeric(object@data)) return("\n data matrix must be numeric")
      if(!all(object@data>=0))return("\n all data matrix elements must be greater or equal to zero")
      if(is.null(rownames(object@data))) return("\n data matrix must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@data))) return("\n data matrix must have the cols names where the taxa IDs were saved")
      if(any(duplicated(rownames(object@data)))) return("\n find in data matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@data)))) return("\n find in data matrix at least a duplicated col name / taxa ID.")
    }
    
    #CHECK META_SAMPLE
    #-------------------------------------#
    if( length(object@meta_sample)!=0 ){
      if(is.null(rownames(object@meta_sample))) return("\n meta_sample data.frame must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@meta_sample))) return("\n meta_sample data.frame have the cols names where the experimental variables were saved")
      if(any(duplicated(rownames(object@meta_sample)))) return("\n find in meta_sample matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@meta_sample)))) return("\n find in meta_sample matrix at least a duplicated col name / experimental variable.")
    }
    
    #CHECK TAXA
    #-------------------------------------#
    if( length(object@taxa)!=0 ){
      if(!is.character(object@taxa)) return("\n taxa matrix must be character")
      if(!all(validUTF8(object@taxa))) return("\n all taxa matrix elements must be encoded with UTF-8")
      if(is.null(rownames(object@taxa))) return("\n taxa matrix must have the rows names where the taxa IDs were saved.")
      if(is.null(colnames(object@taxa))) return("\n taxa matrix must have the cols names where the taxonomic ranks were saved")
      if(any(duplicated(rownames(object@taxa)))) return("\n find in taxa matrix at least a duplicated row name / taxa ID.")
      if(any(duplicated(colnames(object@taxa)))) return("\n find in taxa matrix at least a duplicated col name / rank.")
      if(any(duplicated(rownames(object@taxa[,ncol(object@taxa)])))) return("\n find in last column taxa matrix at least a duplicated taxa ID.")
    }
    
    #CHECK META_TAXA
    #-------------------------------------#
    if( length(object@meta_taxa)!=0 ){
      if(is.null(rownames(object@meta_taxa))) return("\n meta_taxa data.frame must have the rows names where the taxa IDs were saved.")
      if(is.null(colnames(object@meta_taxa))) return("\n meta_taxa data.frame have the cols names where the additional taxa info were saved")
      if(any(duplicated(rownames(object@meta_taxa)))) return("\n find in meta_taxa matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@meta_taxa)))) return("\n find in meta_taxa matrix at least a duplicated col name / additional taxa info.")
    }
    
    # CHECK NETW
    #-------------------------------------#
    if(length(object@netw)!=0){
      if(is_directed(object@netw)) return("\n netw slot must be an undirected igraph object")
      if(!is_weighted(object@netw)) return("\n netw slot must be an weighted igraph object")
      if(is.null(V(g)$name)) return("\n netw vertices name ( V(netw(obj))$name ) cannot be empty.")
    }
    
    # CHECK COMM
    #-------------------------------------#
    if(length(object@comm)!=0 & length(object@netw)==0) return("\n you can't have the comm slot without its associated netw")
    if(length(object@comm)!=0){
      if(is.null(object@comm$name)) return("\n missing vertex names in communities")
      if(object@comm$vcount!=vcount(object@netw)) stop("comm and netw must have the same number of vertices")
      if(!all(object@comm$names==V(object@netw)$name)) stop("communities names must be identical to vertices names in netw")
    }
    
    #CHECK DATA PROPERITES RESPECT OTHER SLOTS
    #-------------------------------------#
    if(length(object@data)!=0){
      
      if(length(object@meta_sample)!=0){
        if(nrow(object@data)!=nrow(object@meta_sample)) return("\n different number of samples in data and meta_sample slots")
        if(!all(rownames(object@data)==rownames(object@meta_sample))) return("\n rows names / sample IDs must be identical in data and meta_sample slots")
      }
      
      if(length(object@taxa)!=0){
        if(ncol(object@data)!=nrow(object@taxa)) return("\n different number of taxa in data and taxa slots")
        if(!all(colnames(object@data)==rownames(object@taxa))) return("\n data colnames and taxa rownames (taxa IDs) must be identical in data and taxa slots")
      }
      
      if(length(object@meta_taxa)!=0){
        if(ncol(object@data)!=nrow(object@meta_taxa)) return("\n different number of taxa in data and meta_taxa slots")
        if(!all(colnames(object@data)==rownames(object@meta_taxa))) return("\n data colnames and meta_taxa rownames (taxa IDs) must be identical in data and taxa slots")
      }
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=ncol(object@data)) return("\n netw vertices number must be equal to columns number of data")
        if(any(colnames(object@data)!=V(object@netw)$name)) return("\n netw vertex name and colnames of data must be equal")
      }
    }
    
    #CHECK TAXA PROPERITES RESPECT OTHER SLOTS
    #-------------------------------------#
    if(length(object@taxa)!=0){
      
      if(length(object@meta_taxa)!=0){
        if(nrow(object@taxa)!=nrow(object@meta_taxa)) return("\n different number of taxa in taxa and meta_taxa slots")
        if(!all(rownames(object@taxa)==rownames(object@meta_taxa))) return("\n taxa rownames and meta_taxa rownames (taxa IDs) must be identical in taxa and taxa slots")
      }
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=nrow(object@taxa)) return("\n netw vertices number must be equal to rows number of taxa")
        if(any(rownames(object@taxa)!=V(object@netw)$name)) return("\n netw vertex name and rownames of taxa must be equal")
      }
    }
    
    #CHECK META_TAXA PROPERITES RESPECT OTHER SLOTS
    #-------------------------------------#
    if(length(object@meta_taxa)!=0){
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=nrow(object@meta_taxa)) return("\n netw vertices number must be equal to rows number of meta_taxa")
        if(any(rownames(object@meta_taxa)!=V(object@netw)$name)) return("\n netw vertex name and rownames of meta_taxa must be equal")
      }
    }
    
    
    TRUE
  })
################################################################################
################################################################################
# END CLASS MG
################################################################################
################################################################################
