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
#' @slot log_data matrix contain the log-ratio transformed data (compostional trans).
#' @slot netw undirected, weighted igraph network.
#' @slot comm graph communities of netw.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy V vcount is_named
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
         log_data="matrix",
         netw="igraph",
         comm="communities"
         ),
  
  prototype=prototype(data=matrix(nrow=0,ncol=0),
                      meta_sample=data.frame(),
                      taxa=matrix(nrow=0,ncol=0),
                      meta_taxa=data.frame(),
                      log_data=matrix(nrow=0, ncol=0),
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
      if(any(is.na(object@taxa))) return("\n in taxa cannot be present NA")
      if(any(object@taxa=="")) return("\n in taxa cannot be present empty character elements (can you rename as unknow or somenthing similar??)")
    }
    
    #CHECK META_TAXA
    #-------------------------------------#
    if( length(object@meta_taxa)!=0 ){
      if(is.null(rownames(object@meta_taxa))) return("\n meta_taxa data.frame must have the rows names where the taxa IDs were saved.")
      if(is.null(colnames(object@meta_taxa))) return("\n meta_taxa data.frame have the cols names where the additional taxa info were saved")
      if(any(duplicated(rownames(object@meta_taxa)))) return("\n find in meta_taxa matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@meta_taxa)))) return("\n find in meta_taxa matrix at least a duplicated col name / additional taxa info.")
    }
    
    #CHECK LOG_DATA
    #-------------------------------------#
    if( length(object@log_data)!=0 ){
      if(!is.numeric(object@log_data)) return("\n log_data matrix must be numeric")
      if(is.null(rownames(object@log_data))) return("\n log_data matrix must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@log_data))) return("\n log_data matrix must have the cols names where the taxa IDs were saved")
      if(any(duplicated(rownames(object@log_data)))) return("\n find in log_data matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@log_data)))) return("\n find in log_data matrix at least a duplicated col name / taxa ID.")
    }
    
    # CHECK NETW
    #-------------------------------------#
    if(length(object@netw)!=0){
      if(!igraph::is_named(object@netw)) return("\n netw vertices must be named")
      if(is_directed(object@netw)) return("\n netw slot must be an undirected igraph object")
      if(!is_weighted(object@netw)) return("\n netw slot must be an weighted igraph object")
      if(any(diag(as_adjacency_matrix(object@netw,sparse=F))!=0)) return("\n netw cannot has self loops")
    }
    
    # CHECK COMM
    #-------------------------------------#
    if(length(object@comm)!=0 & length(object@netw)==0) return("\n you can't have the comm slot without its associated netw")
    if(length(object@comm)!=0){
      if(object@comm$vcount!=vcount(object@netw)) stop("comm and netw must have the same number of vertices")
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
      
      if(length(object@log_data)!=0){
        if(nrow(object@data)!=nrow(object@log_data)) return("\n different number of samples in data and log_data slots")
        if(!all(rownames(object@data)==rownames(object@log_data))) return("\n rows names / sample IDs must be identical in data and log_data slots")
        if(ncol(object@data)!=ncol(object@log_data)) return("\n different number of taxa in data and log_data slots")
        if(!all(colnames(object@data)==colnames(object@log_data))) return("\n cols names / taxa IDs must be identical in data and log_data slots")
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
        if(!all(rownames(object@taxa)==rownames(object@meta_taxa))) return("\n taxa rownames and meta_taxa rownames (taxa IDs) must be identical.")
      }
      
      if(length(object@log_data)!=0){
        if(nrow(object@taxa)!=ncol(object@log_data)) return("\n different number of taxa in taxa and log_data slots")
        if(!all(rownames(object@taxa)==colnames(object@log_data))) return("\n taxa rownames and log_data colnames (taxa IDs) must be identical.")
      }
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=nrow(object@taxa)) return("\n netw vertices number must be equal to rows number of taxa")
        if(any(rownames(object@taxa)!=V(object@netw)$name)) return("\n netw vertex name and rownames of taxa must be equal")
      }
    }
    
    #CHECK META_TAXA PROPERITES RESPECT OTHER SLOTS
    #-------------------------------------#
    if(length(object@meta_taxa)!=0){
      
      if(length(object@log_data)!=0){
        if(nrow(object@meta_taxa)!=ncol(object@log_data)) return("\n different number of taxa in meta_taxa and log_data slots")
        if(!all(rownames(object@meta_taxa)==colnames(object@log_data))) return("\n meta_taxa rownames and log_data colnames (taxa IDs) must be identical.")
      }
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=nrow(object@meta_taxa)) return("\n netw vertices number must be equal to rows number of meta_taxa")
        if(any(rownames(object@meta_taxa)!=V(object@netw)$name)) return("\n netw vertex name and rownames of meta_taxa must be equal")
      }
    }
    
    #CHECK LOG_DATA PROPERITES RESPECT OTHER SLOTS
    #-------------------------------------#
    if(length(object@log_data)!=0){
      
      if(length(object@netw)!=0){
        if(vcount(object@netw)!=ncol(object@log_data)) return("\n netw vertices number must be equal to cols number of log_data")
        if(any(colnames(object@log_data)!=V(object@netw)$name)) return("\n netw vertex name and colnames of log_data must be equal")
      }
    }
    
    TRUE
  })
################################################################################
################################################################################
# END CLASS MG
################################################################################
################################################################################




# CHECK MGNET LIST INTERNAL
################################################################################
#' Check mgnet list
#' 
#' @description Internal function to check if a list is compose from only 
#' mgnet objects.
#' @param object list 
#' @keywords internal
is_mgnet_list <- function(L){
  if(!is.list(L)) stop("the argument must be a list")
  check <- sapply(L, function(x) inherits(x,"mgnet"))
  
  if(any(check==FALSE)){
    stop("list must contain all mgnet objects")
  } else {
    TRUE
  }
}
################################################################################




################################################################################
################################################################################
# SAVE SAMPLE SUM
################################################################################
################################################################################
#' Save sample sum in mgnet.
#' 
#' @description 
#' Store in mgnet object the variable sample_sum in meta_sample in which are
#' saved the count sum for each sample. It is necessary to retrieves relative
#' abundances of samples.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname save_sample_sum-methods
setGeneric("save_sample_sum", function(object) standardGeneric("save_sample_sum"))
#' @rdname save_sample_sum-methods
setMethod("save_sample_sum", "mgnet", function(object){
  
  if(length(object@data)==0) stop("data cannot be empty")
  
  if(length(object@meta_sample)==0){
    object@meta_sample <- data.frame("sample_sum"=rowSums(object@data))
  } else if(length(object@meta_sample)!=0 & !("sample_sum"%in%colnames(object@meta_sample))){
    object@meta_sample$sample_sum <- rowSums(object@data)
  } else {
    object@meta_sample$sample_sum <- rowSums(object@data)
    message("the sample_sum column of meta_sample slot has been rewritten")
  }
  
  return(object)
})
#' @rdname save_sample_sum-methods
setMethod("save_sample_sum","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="save_sample_sum",signature="mgnet"))})
################################################################################
################################################################################
# SAVE SAMPLE SUM
################################################################################
################################################################################




################################################################################
################################################################################
# SAVE LOG-TRANSFORMED DATA
################################################################################
################################################################################
#' Save log-transformed data in mgnet.
#' 
#' @description 
#' Store in mgnet the log-transformed data to dealing with compositional issues.
#'
#' @param object mgnet or list 
#' @param method character with possible values "base" or "inter-quantile".
#' @param zero.dealing character with possible values "plus", "subs" or "none".
#'
#' @export
#' @docType methods
#' @rdname save_log_data-methods
setGeneric("save_log_data", function(object, method="base", zero.dealing="plus") standardGeneric("save_log_data"))
#' @rdname save_log_data-methods
setMethod("save_log_data", "mgnet", function(object, method="base", zero.dealing="plus"){
  
  method <- match.arg(method, c("base","inter-quantile"))
  zero.dealing <- match.arg(zero.dealing, c("none","plus","subs"))
  
  if(zero.dealing!="none"){
    x <- zero_dealing(object@data,mar=1,type=zero.dealing)
  } else {
    x <- object@data
  }
  
  ifelse(method=="base",
         log_x <- clr(x),
         log_x <- iqclr(x))
  
  if(length(log_data(object))!=0) message("log_data slot has been rewritten")
  
  return(mgnet(data=object@data, meta_sample=object@meta_sample,
               taxa=object@taxa, meta_taxa=object@meta_taxa,
               netw=object@netw, comm=object@comm,
               log_data=log_x))
})
#' @rdname save_log_data-methods
setMethod("save_log_data","list",
          function(object, method){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="save_log_data",signature="mgnet"),
                   method=method)})
################################################################################
################################################################################
# SAVE LOG-TRANSFORMED DATA
################################################################################
################################################################################




################################################################################
################################################################################
# CONSTRUCTOR MGNET
################################################################################
################################################################################
#' User constructor for mgnet s4 class.
#' 
#' @description User constructor to create an object belonging to formal s4 
#' class mgnet avoiding the new function.
#' 
#' @param data numeric matrix with all elements >=0.  
#' @param meta_sample data.frame with experimental variables.
#' @param taxa character matrix with taxonomic classification.  
#' @param meta_taxa data.frame with addition info on taxa.
#' @param log_data log-transformed data.
#' @param netw un-directed weighted igraph network or the associated adjacency matrix.
#' @param comm communities object associated to netw.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
mgnet <- function(data=matrix(nrow=0,ncol=0),
                  meta_sample=data.frame(),
                  taxa=matrix(nrow=0,ncol=0),
                  meta_taxa=data.frame(),
                  log_data=matrix(nrow=0,ncol=0),
                  netw=make_empty_graph(n=0, directed=FALSE),
                  comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))
){
  
  if(length(netw)!=0){
    if(is.matrix(netw)){
      if(!is.numeric(netw) | !is.matrix(netw) | !isSymmetric(netw)) stop("netw if it is given as adjacency matrix must be simmetric")
      netw <- igraph::graph_from_adjacency_matrix(netw,'undirected',weighted=TRUE)
    }
  }
  
  if(length(data)!=0 & length(meta_sample)==0){
    meta_sample <- data.frame("sample_sum"=rowSums(data))
  } else if( length(meta_sample)!=0 & !("sample_sum"%in%colnames(meta_sample))){
    meta_sample$sample_sum <- rowSums(data)
  }
  
  return(new("mgnet",data=data, meta_sample=meta_sample,
             taxa=taxa, meta_taxa=meta_taxa,
             log_data=log_data,
             netw=netw, comm=comm))
  
}
################################################################################
################################################################################
# END CONSTRUCTOR MGNET
################################################################################
################################################################################






################################################################################
################################################################################
# GETTERS MG
################################################################################
################################################################################
# DATA
#####################################
#' Retrieves ngs data.
#' 
#' @description 
#' Return the numeric matrix associated with the ngs data from mgnet class.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname data-methods
setGeneric("data", function(object) standardGeneric("data"))
#' @rdname data-methods
setMethod("data", "mgnet", function(object) object@data)
#' @rdname data-methods
setMethod("data","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="data",signature="mgnet"))})
#####################################
# META SAMPLE
#####################################
#' Retrieves sample metadata.
#'
#' @description
#' Return the data.frame associated with the sample metadata with experimental
#' variables from mg class.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname meta_sample-methods
setGeneric("meta_sample", function(object) standardGeneric("meta_sample"))
#' @rdname meta_sample-methods
setMethod("meta_sample", "mgnet", function(object) object@meta_sample)
#' @rdname meta_sample-methods
setMethod("meta_sample","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="meta_sample",signature="mgnet"))})
#####################################
# TAXA
#####################################
#' Retrieves taxonomy table.
#'
#' @description
#' Return the matrix associated with taxonomy classification from mg class.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname taxa-methods
setGeneric("taxa", function(object) standardGeneric("taxa"))
#' @rdname taxa-methods
setMethod("taxa", "mgnet", function(object) object@taxa)
#' @rdname taxa-methods
setMethod("taxa","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="taxa",signature="mgnet"))})
#####################################
# META TAXA
#####################################
#' Retrieves taxa additional info.
#'
#' @description
#' Return the data.frame associated with the additional taxa information from 
#' mg class.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname meta_taxa-methods
setGeneric("meta_taxa", function(object) standardGeneric("meta_taxa"))
#' @rdname meta_taxa-methods
setMethod("meta_taxa", "mgnet", function(object) object@meta_taxa)
#' @rdname meta_taxa-methods
setMethod("meta_taxa","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="meta_taxa",signature="mgnet"))})
#####################################
# LOG DATA
#####################################
#' Retrieves log-transformed data.
#'
#' @description
#' Return the matrix associated with the log-transformed data 
#' mg class.
#'
#' @param object mgnet or list 
#'
#' @export
#' @docType methods
#' @rdname log_data-methods
setGeneric("log_data", function(object) standardGeneric("log_data"))
#' @rdname log_data-methods
setMethod("log_data", "mgnet", function(object) object@log_data)
#' @rdname log_data-methods
setMethod("log_data","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="log_data",signature="mgnet"))})
#####################################
# NETW
#####################################
#' Retrieves netw.
#' 
#' @description 
#' Return the igraph un-directed weighted network.
#'
#' @usage netw(object)
#'
#' @param object mgnet or list.
#'
#' @export
#' @docType methods
#' @rdname netw-methods
setGeneric("netw", function(object) standardGeneric("netw"))
#' @rdname netw-methods
setMethod("netw", "mgnet", function(object){return(object@netw)})
#' @rdname netw-methods
setMethod("netw","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="netw",signature="mgnet"))})
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
#' @param object mgnet or list.
#'
#' @export
#' @docType methods
#' @rdname comm-methods
setGeneric("comm", function(object) standardGeneric("comm"))
#' @rdname comm-methods
setMethod("comm", "mgnet", function(object){return(object@comm)})
#' @rdname comm-methods
setMethod("comm","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="comm",signature="mgnet"))})
#####################################
# ADJACENCY MATRIX
#####################################
#' Retrieves adjacency matrix.
#'
#' @description
#' Return the adjacency matrix associated to netw slot using the igraph function 
#' as_adjacency_matrix.
#'
#' @usage adjacency_matrix(object)
#'
#' @param object mgnet or list.
#'
#' @importFrom igraph as_adjacency_matrix
#' @export
#' @docType methods
#' @rdname adjacency_matrix-methods
setGeneric("adjacency_matrix", function(object) standardGeneric("adjacency_matrix"))
#' @rdname adjacency_matrix-methods
setMethod("adjacency_matrix", "mgnet", function(object){
  if(length(netw(object))!=0){
    return(igraph::as_adjacency_matrix(netw(object),attr="weight",sparse=F))
    } else {
    return(matrix(nrow=0,ncol=0))
  }})
#' @rdname adjacency_matrix-methods
setMethod("adjacency_matrix","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="adjacency_matrix",signature="mgnet"))})
#####################################
# ADJACENCY LIST
#####################################
#' Retrieves adjacency edge list
#'
#' @description
#' Return the adjacency list associated to netw slot as data.frame using the 
#' igraph function as_data_frame.
#'
#' @usage adjacency_list(object)
#'
#' @param object mgnet or list.
#'
#' @importFrom igraph as_data_frame
#' @export
#' @docType methods
#' @rdname adjacency_list-methods
setGeneric("adjacency_list", function(object) standardGeneric("adjacency_list"))
#' @rdname adjacency_list-methods
setMethod("adjacency_list", "mgnet", function(object){
  if(length(netw(object))!=0){
    return(igraph::as_data_frame(netw(object),what="edges"))
  } else {
    stop("there isn't netw")
  }})
#' @rdname adjacency_list-methods
setMethod("adjacency_list","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="adjacency_list",signature="mgnet"))})
################################################################################
################################################################################
# END GETTERS MGNET
################################################################################
################################################################################




# CHECK ASSIGN LISTS
################################################################################
#' Check lists in assign methods
#' 
#' @description Internal function to check if the two list in assign methods are
#' coherent.
#' @param object mgnet list
#' @param value list 
#' @keywords internal
are_lists_assign <- function(object, value){
  if(!is.list(object) | !is.list(value)) stop("object and value must be a lists")
  if(length(object)!=length(value)) stop("lengths of object and value lists must be equal")
  if(is.null(names(object)) | is.null(names(value))) stop("object and value list have not named element")
  if(any(names(object)!=names(value))) stop("object and value elements must have the same named element in the same order")
  
  # Check object is a list of mgnet
  check <- sapply(object, function(x) inherits(x,"mgnet"))
  if(any(check==FALSE)){
    stop("object must contain all mgnet objects")
  } else {
    TRUE
  }
}
################################################################################




################################################################################
################################################################################
# SETTERS MGNET
################################################################################
################################################################################
# DATA<-
#####################################
#' Assign a new ngs matrix to \code{object}
#'
#' @usage data(object) <- value
#'
#' @param object mgnet.
#' @param value matrix
#'
#' @docType methods
#' @rdname assign-data
#' @export
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
#' @rdname assign-data
setMethod("data<-", c("mgnet", "matrix"), function(object, value){
  new("mgnet",data=value, meta_sample=object@meta_sample, 
      taxa=object@taxa, meta_taxa=object@meta_taxa,
      log_data=object@log_data,
      netw=object@netw, comm=object@comm)
})
#####################################
# META SAMPLE<-
#####################################
#' Assign new samples metadata to \code{object}
#'
#' @usage meta_sample(object) <- value
#'
#' @param object mgnet
#' @param value data.frame
#'
#' @export
#' @docType methods
#' @rdname assign-meta_sample
setGeneric("meta_sample<-", function(object, value) standardGeneric("meta_sample<-"))
#' @rdname assign-meta_sample
setMethod("meta_sample<-", c("mgnet", "data.frame"), function(object, value){
  new("mgnet",data=object@data, meta_sample=value, 
      taxa=object@taxa, meta_taxa=object@meta_taxa,
      log_data=object@log_data,
      netw=object@netw, comm=object@comm)
})
#####################################
# TAXA<-
#####################################
#' Assign new taxonomy table to \code{object}
#'
#' @usage taxa(object) <- value
#'
#' @param object mgnet
#' @param value matrix
#'
#' @export
#' @docType methods
#' @rdname assign-taxa
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))
#' @rdname assign-taxa
setMethod("taxa<-", c("mgnet", "matrix"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=value, meta_taxa=object@meta_taxa,
      log_data=object@log_data,
      netw=object@netw, comm=object@comm)
})
#####################################
# META TAXA<-
#####################################
#' Assign new meta taxa to \code{object}
#'
#' @usage meta_taxa(object) <- value
#'
#' @param object mgnet
#' @param value data.frame
#'
#' @export
#' @docType methods
#' @rdname assign-meta_taxa
setGeneric("meta_taxa<-", function(object, value) standardGeneric("meta_taxa<-"))
#' @rdname assign-meta_taxa
setMethod("meta_taxa<-", c("mgnet", "data.frame"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=object@taxa, meta_taxa=value,
      log_data=object@log_data,
      netw=object@netw, comm=object@comm)
})
#####################################
# LOG DATA<-
#####################################
#' Assign new log_data to \code{object}
#'
#' @usage log_data(object) <- value
#'
#' @param object mgnet
#' @param value matrix
#'
#' @export
#' @docType methods
#' @rdname assign-log_data
setGeneric("log_data<-", function(object, value) standardGeneric("log_data<-"))
#' @rdname assign-log_data
setMethod("log_data<-", c("mgnet", "matrix"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=object@taxa, meta_taxa=object@meta_taxa,
      log_data=value,
      netw=object@netw, comm=object@comm)
})
#####################################
# NETW<-
#####################################
#' Assign a new netw to \code{object}
#'
#' @usage netw(object) <- value
#'
#' @param object mgnet.
#' @param value \code{\link{igraph}} or the equivalent adjacency matrix.
#'
#' @export
#' @docType methods
#' @rdname assign-netw
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))
#' @rdname assign-netw
setMethod("netw<-", "mgnet", function(object, value){

  if( !is.igraph(value) & !is.matrix(value) ) stop("value must be a matrix or an igraph object")
  
  mgnet(data=object@data, meta_sample=object@meta_sample, 
        taxa=object@taxa, meta_taxa=object@meta_taxa,
        log_data=object@log_data,
        netw=value, comm=object@comm)
})
#####################################
# ADJACENCY LIST<-
#####################################
#' Assign a new netw to \code{object} from an adjacency edge data.frame
#'
#' @usage adjacency_list(object) <- value
#'
#' @param object mgnet.
#' @param value \code{\link{data.frame}} with three named columns "from","to","weight"
#' like in the as_data_frame function of igraph.
#'
#' @importFrom igraph make_clusters cluster_fast_greedy graph_from_adjacency_matrix add_edges
#'
#' @export
#' @docType methods
#' @rdname assign-adjacency_list
setGeneric("adjacency_list<-", function(object, value) standardGeneric("adjacency_list<-"))
#' @rdname assign-adjacency_list
setMethod("adjacency_list<-", c("mgnet", "data.frame"), function(object, value){
  
  if(!is.data.frame(value)) stop("value must be a data.frame")
  if(!is.numeric(value[,3])) stop("value column must be numeric")
  if(!any(value[,1]%in%taxaID(object))) stop("find at leat a vertex in value not present in object")
  if(!any(value[,2]%in%taxaID(object))) stop("find at leat a vertex in value not present in object")
  if(ncol(value)!=3) stop("value must be a data.frame with three named columns [from,to,weight]")
  if(any(colnames(value)!=c("from","to","weight"))) stop("value must be a data.frame with three named columns [from,to,weight]")
  
  empty_adj <- matrix(0,nrow=ntaxa(object), ncol=ntaxa(object),
                      dimnames=list(taxaID(object),taxaID(object)))
  new_netw <- igraph::graph_from_adjacency_matrix(empty_adj, mode="undirected",weighted=T)
  
  for(e in 1:nrow(value)){
    new_netw <- igraph::add_edges(new_netw, c(value[e,1],value[e,2]), weight=value[e,3])
  }
  
  if(length(netw(object))!=0) message("netw slot has been rewritten")
  mgnet(data=object@data, meta_sample=object@meta_sample, 
        taxa=object@taxa, meta_taxa=object@meta_taxa,
        log_data=object@log_data,
        netw=new_netw,
        comm=igraph::cluster_fast_greedy(igraph::make_empty_graph(n=0, directed=FALSE)))
})
#' @rdname assign-adjacency_list
setMethod("adjacency_list<-",c("list","list"),
          function(object,value){
            are_lists_assign(object,value)
            for(i in 1:length(object)){
              adjacency_list(object[[i]]) <- value[[i]]}
            return(object)})
#####################################
# COMM<-
#####################################
#' Assign a new comm to \code{object}
#'
#' @usage comm(object) <- value
#'
#' @param object mgnet.
#' @param value communities
#' 
#' @importFrom igraph make_clusters 
#' 
#' @export
#' @docType methods
#' @rdname assign-comm
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))
#' @rdname assign-comm
setMethod("comm<-", c("mgnet", "communities"), function(object, value){
  
  if(is.null(names(value$membership))){
    
    if(length(object@netw)==0) stop("netw slot in object cannot be empty")
    if(ntaxa(object)!=length(value$membership)) stop("vertices number in object and in value is different")
    
    return(mgnet(data=object@data, meta_sample=object@meta_sample,
                 taxa=object@taxa, meta_taxa=object@meta_taxa,
                 log_data=object@log_data,
                 netw=object@netw, comm=value))
    
  } else {
    
    if(length(netw(object))==0) stop("netw slot cannot be empty")
    if(all(taxaID(object)%in%names(value$membership))) stop("all taxaID in value must be present also in object")
    
    namedMemb <- value$membership
    namedMemb.new <- setNames(rep(0,ntaxa(object)), taxaID(object))
    namedMemb.new[names(namedMemb)] <- namedMemb
    
    comm.new <- igraph::make_clusters(object@netw,
                                      membership=as.numeric(namedMemb.new),
                                      algorithm="NA",
                                      modularity=1)
    comm.new$modularity <- NA
    
    return(mgnet(data=object@data, meta_sample=object@meta_sample,
                 taxa=object@taxa, meta_taxa=object@meta_taxa,
                 log_data=object@log_data,
                 netw=object@netw, comm=comm.new))
    
  }
})
#' @rdname assign-comm
setMethod("comm<-",c("list","list"),
          function(object,value){
            are_lists_assign(object,value)
            classes <- sapply(value,class)
            
            if( any(classes!="communities") ){
              stop("all elements in value list must belong to communities or mgnet class")
            }
            
            for(i in 1:length(object)){
              comm(object[[i]]) <- value[[i]]}
            return(object)})
################################################################################
################################################################################
# END SETTERS MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# BASE METHODS
################################################################################
################################################################################
# NSAMPLE
#####################################
#' Get number of samples.
#' 
#' @description 
#' Return an integer indicating the number of sample.
#'
#' @usage nsample(object)
#'
#' @param object mgnet ot list.
#'
#' @rdname nsample-methods
#' @docType methods
#' @export
setGeneric("nsample", function(object) standardGeneric("nsample"))
#' @rdname nsample-methods
setMethod("nsample", "mgnet", function(object){
  if(length(object@data!=0)) return(nrow(object@data))
  else if(length(object@meta_sample!=0)) return(nrow(object@meta_sample))
  else if(length(object@log_data!=0)) return(nrow(object@log_data))
  else return(NULL)
})
#' @rdname nsample-methods
setMethod("nsample","list",
          function(object){
            is_mgnet_list(object)
            sapply(object, selectMethod(f="nsample",signature="mgnet"))})
#####################################
# NTAXA
#####################################
#' Get number of taxa
#' 
#' @description 
#' Return an integer indicating the number of taxa
#'
#' @usage ntaxa(object)
#'
#' @param object mgnet or list.
#'
#' @importFrom igraph vcount
#' @rdname ntaxa-methods
#' @docType methods
#' @export
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))
#' @rdname ntaxa-methods
setMethod("ntaxa", "mgnet", function(object){
  if(length(object@data)!=0) return(ncol(object@data))
  else if(length(object@taxa)!=0) return(nrow(object@taxa))
  else if(length(object@meta_taxa)!=0) return(nrow(object@meta_taxa))
  else if(length(object@log_data)!=0) return(ncol(object@log_data))
  else if(length(object@netw)!=0) return(vcount(object@netw))
  else return(NULL)
})
#' @rdname ntaxa-methods
setMethod("ntaxa","list",
          function(object){
            is_mgnet_list(object)
            sapply(object, selectMethod(f="ntaxa",signature="mgnet"))})
#####################################
# SAMPLE NAME
#####################################
#' Get samples names.
#' 
#' @description 
#' Return names of samples as character vector.
#'
#' @usage sample_name(object)
#'
#' @param object mgnet or list.
#'
#' @rdname sample_name-methods
#' @docType methods
#' @export
setGeneric("sample_name", function(object) standardGeneric("sample_name"))
#' @rdname sample_name-methods
setMethod("sample_name", "mgnet", function(object){
  if(length(object@data!=0)) return(rownames(object@data))
  else if(length(object@meta_sample!=0)) return(rownames(object@meta_sample))
  else if(length(object@log_data!=0)) return(rownames(object@log_data))
  else return(NULL)
})
#' @rdname sample_name-methods
setMethod("sample_name","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="sample_name",signature="mgnet"))})
#####################################
# TAXA ID
#####################################
#' Get taxa ID.
#' 
#' @description 
#' Return taxonomy id as a character vector.
#'
#' @usage taxaID(object)
#'
#' @param object mgnet or list.
#'
#' @importFrom igraph V
#' @rdname taxaID-methods
#' @docType methods
#' @export
setGeneric("taxaID", function(object) standardGeneric("taxaID"))
#' @rdname taxaID-methods
setMethod("taxaID", "mgnet", function(object){
  if(length(object@data!=0)) return(colnames(object@data))
  else if(length(object@taxa!=0)) return(rownames(object@taxa))
  else if(length(object@meta_taxa)!=0) return(rownames(object@meta_taxa))
  else if(length(object@log_data)!=0) return(colnames(object@log_data))
  else if(length(object@netw)!=0) return(V(object@netw)$name)
  else return(0)
})
#' @rdname taxaID-methods
setMethod("taxaID","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="taxaID",signature="mgnet"))})
#####################################
# RANKS 
#####################################
#' Get taxonomic ranks.
#' 
#' @description 
#' Return taxonomy ranks as a character vector.
#'
#' @usage ranks(object)
#'
#' @param object mgnet or list.
#'
#' @rdname ranks-methods
#' @docType methods
#' @export
setGeneric("ranks", function(object) standardGeneric("ranks"))
#' @rdname ranks-methods
setMethod("ranks", "mgnet", function(object){return(colnames(object@taxa))})
#' @rdname ranks-methods
setMethod("ranks","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="ranks",signature="mgnet"))})
#####################################
# NRANKS 
#####################################
#' Get taxonomic ranks number.
#' 
#' @description 
#' Return taxonomy ranks number number as integer.
#'
#' @usage nrank(object)
#'
#' @param object mgnet or list.
#'
#' @rdname nrank-methods
#' @docType methods
#' @export
setGeneric("nrank", function(object) standardGeneric("nrank"))
#' @rdname nrank-methods
setMethod("nrank", "mgnet", function(object){return(ncol(object@taxa))})
#' @rdname nrank-methods
setMethod("nrank","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="nrank",signature="mgnet"))})
#####################################
# SAMPLE_INFO 
#####################################
#' Get sample metadata variables.
#' 
#' @description 
#' Return sample metadata variables number as character vector.
#'
#' @usage sample_info(object)
#'
#' @param object mgnet or list.
#'
#' @rdname sample_info-methods
#' @docType methods
#' @export
setGeneric("sample_info", function(object) standardGeneric("sample_info"))
#' @rdname sample_info-methods
setMethod("sample_info", "mgnet", function(object){return(colnames(object@meta_sample))})
#' @rdname sample_info-methods
setMethod("sample_info","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="sample_info",signature="mgnet"))})
#####################################
# TAXA_INFO
#####################################
#' Get taxa additional info.
#' 
#' @description 
#' Return taxa additional info as character vector.
#'
#' @usage taxa_info(object)
#'
#' @param object mgnet or list.
#'
#' @rdname taxa_info-methods
#' @docType methods
#' @export
setGeneric("taxa_info", function(object) standardGeneric("taxa_info"))
#' @rdname taxa_info-methods
setMethod("taxa_info", "mgnet", function(object){return(colnames(object@meta_taxa))})
#' @rdname taxa_info-methods
setMethod("taxa_info","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="taxa_info",signature="mgnet"))})
#####################################
# TAXA NAME 
#####################################
#' Get taxa name.
#' 
#' @description 
#' Get taxa name at choosen rank as character vector.
#' 
#' @usage taxa_name(object, rank)
#' 
#' @param object mgnet or list.
#' @param rank taxonomic level choosen (if it is not set, the finest taxonomic rank is assumed)
#' 
#' @rdname taxa_name-methods
#' @docType methods
#' @export
setGeneric("taxa_name", function(object, rank) standardGeneric("taxa_name"))
#' @rdname taxa_name-methods
setMethod("taxa_name", c("mgnet","missing"), function(object) object@taxa[,nrank(object)])
#' @rdname taxa_name-methods
setMethod("taxa_name", c("mgnet","character"),
          function(object, rank){
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",
                                                    toString(ranks(object)),"}"))
            return(object@taxa[,rank])
          })
#' @rdname taxa_name-methods
setMethod("taxa_name",c("list","missing"),
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="taxa_name",signature=c("mgnet","missing")))})
#' @rdname taxa_name-methods
setMethod("taxa_name",c("list","character"),
          function(object,rank){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="taxa_name",signature=c("mgnet","character")),
                   rank=rank)}
)
#####################################
# SAMPLE SUM 
#####################################
#' Get samples sum.
#' 
#' @description 
#' Retrieves the sum on each sample using the sample_sum variable in in 
#' meta_sample slot if it's present (see \code{\link{save_sample_sum}}).
#' 
#' @usage sample_sum(object)
#' 
#' @param object mgnet or list.
#' 
#' @rdname sample_sum-methods
#' @docType methods
#' @export
setGeneric("sample_sum", function(object) standardGeneric("sample_sum"))
#' @rdname sample_sum-methods
setMethod("sample_sum", "mgnet",function(object){
  if(length(object@meta_sample)!=0 & "sample_sum"%in%colnames(object@meta_sample)){
    return(object@meta_sample$sample_sum)
  } else {
    return(NULL)
  }
})
#' @rdname sample_sum-methods
setMethod("sample_sum","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="sample_sum",signature="mgnet"))})
#####################################
# ABUNDANCE 
#####################################
#' Get abundances at choosen rank.
#' 
#' @description 
#' Retrieves the abundances of data at choosen taxonomy rank. Abundance at 
#' higher rank are returned as sums of abundance of elements with the same
#' classification.
#' 
#' @usage abundance(object,rank)
#' 
#' @param object mgnet or list.
#' @param rank (Optional) character with the taxonomic rank choosen.
#' 
#' @importFrom stats aggregate
#' 
#' @export
#' @docType methods
#' @rdname abundance-methods
setGeneric("abundance", function(object,rank) standardGeneric("abundance"))
#' @rdname abundance-methods
setMethod("abundance", c("mgnet","missing"),function(object)return(object@data))
#' @rdname abundance-methods
setMethod("abundance", c("mgnet","character"),function(object,rank){
  
  if(length(object@data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
  if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
  
  different.taxa <- unique(object@taxa[,rank])
  data.aggregate <- t(aggregate(t(object@data), by=list(object@taxa[,rank]), 
                                FUN="sum", drop=FALSE, simplify = TRUE))
  colnames(data.aggregate) <- data.aggregate[1,]
  data.aggregate <- data.aggregate[-1,]
  data.aggregate <- data.aggregate[,different.taxa]
  class(data.aggregate) <- "numeric"
  
  return(data.aggregate)
})
#' @rdname abundance-methods
setMethod("abundance",c("list","missing"),
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="abundance",signature=c("mgnet","missing")))})
#' @rdname abundance-methods
setMethod("abundance",c("list","character"),
          function(object,rank){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="abundance",signature=c("mgnet","character")),
                   rank=rank)}
)
#####################################
# RELVATIVE 
#####################################
#' Get relative abundances.
#' 
#' @description 
#' Retrieves the relative abundances of data, normalized by their sample sum at
#' the taxonomic rank choosen. In mgnet object must be save the variable sample_sum
#' in meta_sample slot (see \code{\link{save_sample_sum}}).
#' 
#' @usage relative(object,rank)
#' 
#' @param object mgnet or mgnet-list.
#' @param rank (Optional) character with the taxonomic rank choosen.
#' 
#' @importFrom stats aggregate
#' 
#' @rdname relative-methods
#' @docType methods
#' @export
setGeneric("relative", function(object,rank) standardGeneric("relative"))
#' @rdname relative-methods
setMethod("relative", c("mgnet","missing"),function(object){
  
  if(length(object@data)==0) stop("data cannot be empty")
  return(object@data/object@meta_sample$sample_sum)
  
})
#' @rdname relative-methods
setMethod("relative", c("mgnet","character"),function(object,rank){
  if(length(object@data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
  if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
  if(!("sample_sum"%in%sample_info(object))) stop("to relative must be present sample_sum in meta_sample (See save_sample_sum).")
  
  different.taxa <- unique(object@taxa[,rank])
  data.aggregate <- t(aggregate(t(object@data), by=list(object@taxa[,rank]), 
                                FUN="sum", drop=FALSE))
  colnames(data.aggregate) <- data.aggregate[1,]
  data.aggregate <- data.aggregate[-1,]
  data.aggregate <- data.aggregate[,different.taxa]
  class(data.aggregate) <- "numeric"
  
  return(data.aggregate/object@meta_sample$sample_sum)
})
#' @rdname relative-methods
setMethod("relative",c("list","missing"),
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="relative",signature=c("mgnet","missing")))})
#' @rdname relative-methods
setMethod("relative",c("list","character"),
          function(object,rank){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="relative",signature=c("mgnet","character")),
                   rank=rank)}
)
#####################################
# NCOMM
#####################################
#' Get communities number.
#' 
#' @description 
#' Return an integer indicates the communities present in comm slot.
#'
#' @usage ncomm(object)
#'
#' @param object mgnet or mgnet-list.
#'
#' @importFrom igraph sizes
#' @rdname ncomm-methods
#' @docType methods
#' @export
setGeneric("ncomm", function(object) standardGeneric("ncomm"))
#' @rdname ncomm-methods
setMethod("ncomm", "mgnet", function(object){
  if(length(object@comm)!=0){
    sizes <- names(igraph::sizes(object@comm))
    class(sizes) <- "numeric"
    return(max(sizes))
  } else {
    return(NULL)
  }
})
#' @rdname ncomm-methods
setMethod("ncomm","list",
          function(object){
            is_mgnet_list(object)
            sapply(object, selectMethod(f="ncomm",signature="mgnet"))})
#####################################
# COMMID
#####################################
#' Get communities ID.
#' 
#' @description 
#' Return vector of character with the community ID for each vertex.
#'
#' @usage commID(object)
#'
#' @param object mgnet or mgnet-list.
#'
#' @importFrom igraph membership
#' @importFrom stats setNames
#' @rdname commID-methods
#' @docType methods
#' @export
setGeneric("commID", function(object) standardGeneric("commID"))
#' @rdname commID-methods
setMethod("commID", "mgnet", function(object){
  if(length(object@comm)!=0){
    return(setNames(as.character(membership(object@comm)),taxaID(object)))
  } else {
    return(NULL)
  }
})
#' @rdname commID-methods
setMethod("commID","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="commID",signature="mgnet"))})
################################################################################
################################################################################
# END BASE METHODS
################################################################################
################################################################################




################################################################################
################################################################################
# SAMPLE SELECTION
################################################################################
################################################################################
#' Manage the communities membership of taxa 
#' 
#' @description The function takes as input a vector of logical or position 
#' indices to evaluate which samples to keep.
#' 
#' @usage selection_sample(object,...,condition="AND")
#' 
#' @param object mgnet.
#' @param ... arbitrary set of vectors of functions. If vectors they must contain
#' numeric or name or logical indices to be preserved. If functions they must return
#' indices as just explained.
#' @param condition (Optional default "AND") logical operation between indices.
#' 
#' @rdname selection_sample-methods
#' @docType methods
#' @export
setGeneric("selection_sample", function(object, ..., condition="AND") standardGeneric("selection_sample"))
#' @rdname selection_sample-methods
setMethod("selection_sample", "mgnet",
          function(object, ..., condition="AND"){
            
            condition <- match.arg(condition, c("AND","OR"))
            if(length(object@netw)!=0) warning("Sample subsetting is not defined for network, the resulting object will not have the slots netw and comm")
            
            IDX <- list(...)
            for(i in 1:length(IDX)){
              
              idx <- IDX[[i]]
              if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
              
              if(is.function(idx)){
                idx <- idx(object)
                if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
              }  
              
              
              if(is.numeric(idx)){
                #numeric
                if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>nsample(object)){
                  stop("numeric indices must be integers in range 1 to the maximum nuber of sample")
                }
                #transform to logical
                IDX[[i]] <- (1:nsample(object))%in%idx
                
              }else if(is.character(idx)){
                #character
                if( !any(idx%in%sample_name(object)) ){
                  stop("string indices must be a subset of sample_name of object")
                }
                #transform to logical
                IDX[[i]] <- sample_name(object)%in%idx
                
              }else if(is.logical(idx)){
                #logical
                if(length(idx)!=nsample(object)){
                  stop("logical indices must have the length equal to sample number")
                }
                IDX[[i]] <- idx
                
              } else {stop("the indices must character, numeric or logical vectors")}
            }
            
            
            IDX <- do.call(rbind,IDX)
            # Condition  
            if(condition=="AND"){
              IDX <- as.logical(apply(IDX,2,prod))
            } else {
              IDX <- apply(IDX,2,sum)>0
            }
            
            ifelse(length(object@data)!=0, data.new<-object@data[IDX,,drop=F], data.new<-object@data)
            ifelse(length(object@meta_sample)!=0, meta.new<-object@meta_sample[IDX,,drop=F], meta.new<-object@meta_sample)
            ifelse(length(object@log_data)!=0, log_data.new<-object@log_data[IDX,,drop=F], log_data.new<-object@log_data)
            
            
            return(mgnet(data=data.new, meta_sample=meta.new,
                         taxa=object@taxa, meta_taxa=object@meta_taxa,
                         log_data=log_data.new))
          })
#' @rdname selection_sample-methods
setMethod("selection_sample","list",
          function(object, ..., condition="AND"){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="selection_sample",signature="mgnet"),
                   ...=..., condition=condition)}
)
################################################################################
################################################################################
# END SAMPLE SELECTION
################################################################################
################################################################################



################################################################################
################################################################################
# TAXA SELECTION
################################################################################
################################################################################
#' Selection a subset of taxa from mgnet object.
#' 
#' @description The function takes as input a vector of logical or position 
#' indices to evaluate which taxa to keep.
#' 
#' @usage selection_taxa(object, ..., condition="AND",trim=TRUE)
#' 
#' @param object mgnet.
#' @param ... arbitrary set of vectors of functions. If vectors they must contain
#' numeric or name or logical indices to be preserved. If functions they must return
#' indices as just explained.
#' @param condition (Optional default "AND") character indicates the logical operation between indices.
#' @param trim (Option default TRUE) logical, if true the discarded taxa are removed from obj otherwise they are only set to zero in data.
#' 
#' @importFrom igraph subgraph
#' @rdname selection_taxa-methods
#' @docType methods
#' @export
setGeneric("selection_taxa", function(object,...,condition="AND",trim=TRUE) standardGeneric("selection_taxa"))
#' @rdname selection_taxa-methods
setMethod("selection_taxa", "mgnet",
          function(object,...,condition="AND",trim=TRUE){
            
            condition <- match.arg(condition, c("AND","OR"))
            if(!is.logical(trim)) stop("trim must be logical")
            
            IDX <- list(...)
            for(i in 1:length(IDX)){
              
              idx <- IDX[[i]]
              if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
              
              if(is.function(idx)){
                idx <- idx(object)
                if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
              }  
              
              if(is.numeric(idx)){
                #numeric
                if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>ntaxa(object)){
                  stop("numeric indices must be integers in range 1 to the maximum nuber of sample")
                }
                #transform to logical
                IDX[[i]] <- (1:ntaxa(object))%in%idx
                
              }else if(is.character(idx)){
                #character
                if( !any(idx%in%taxaID(object)) ){
                  stop("string indices must be a subset of sample_name of object")
                }
                #transform to logical
                IDX[[i]] <- taxaID(object)%in%idx
                
              }else if(is.logical(idx)){
                #logical
                if(length(idx)!=ntaxa(object)){
                  stop("logical indices must have the length equal to sample number")
                }
                IDX[[i]] <- idx
                
              } else {stop("the indices must character, numeric or logical vectors")}
            }
            
            
            IDX <- do.call(rbind,IDX)
            # Condition  
            if(condition=="AND"){
              IDX <- as.logical(apply(IDX,2,prod))
            } else {
              IDX <- apply(IDX,2,sum)>0
            }
            
            if(trim){
              
              ifelse(length(object@data)!=0, data.new<-object@data[,IDX,drop=F], data.new<-object@data)
              ifelse(length(object@taxa)!=0, taxa.new<-object@taxa[IDX,,drop=F], taxa.new<-object@taxa)
              ifelse(length(object@meta_taxa)!=0, meta_taxa.new<-object@meta_taxa[IDX,,drop=F], meta_taxa.new<-object@meta_taxa)
              ifelse(length(object@log_data)!=0, log_data.new<-object@log_data[,IDX,drop=F], log_data.new<-object@log_data)
              
              if(length(object@netw)!=0){
                netw.new<-igraph::subgraph(object@netw,IDX)
              } else {
                netw.new<-object@netw
              }
              
              if(length(object@comm)!=0){
                comm.new <- object@comm
                if(is.character(IDX)) IDX <- which(taxaID(object)%in%IDX)
                comm.new$membership <- object@comm$membership[IDX]
                comm.new$vcount <- length(comm.new$membership)
                comm.new$modularity <- NA
              } else {
                comm.new <- object@comm
              }
              
              return(mgnet(data=data.new,
                           meta_sample=object@meta_sample,
                           taxa=taxa.new,
                           meta_taxa=meta_taxa.new,
                           log_data=log_data.new,
                           netw=netw.new,
                           comm=comm.new))
              
            } else {
              
              if(length(object@data)!=0){
                data.new <- object@data
                data.new[,!IDX] <- 0
              } else {
                data.new <- object@data
              }
              
              if(length(object@log_data)!=0){
                log_data.new <- object@log_data
                log_data.new[,!IDX] <- 0
              } else {
                data.new <- object@log_data
              }
              
              if(length(object@netw)!=0){
                sub.netw <- subgraph(netw(object), taxaID(object)[IDX])
                preserved.edges <- E(object@netw)%in%E(sub.netw)
                netw.new <- subgraph.edges(graph=netw(object), 
                                           eids=E(netw(object))[preserved.edges],
                                           delete.vertices = FALSE)
              } else {
                netw.new<-object@netw
              }
              
              if(length(object@comm)!=0){
                comm.new <- object@comm
                comm.new$membership[!IDX] <- 0
                comm.new$modularity <- NA
              } else {
                comm.new <- object@comm
              }
              
              return(mgnet(data=data.new,
                           meta_sample=object@meta_sample,
                           taxa=object@taxa,
                           meta_taxa=object@meta_taxa,
                           log_data=log_data.new,
                           netw=netw.new,
                           comm=comm.new))
              
            }
            
          })
#' @rdname selection_taxa-methods
setMethod("selection_taxa","list",
          function(object,...,condition="AND",trim=TRUE){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="selection_taxa",signature="mgnet"),
                   ...=..., condition=condition, trim=trim)}
)
################################################################################
################################################################################
# END TAXA SELECTION
################################################################################
################################################################################




################################################################################
################################################################################
# EXTRACTOR MG
################################################################################
################################################################################
#' Method extensions to extraction operator for mg object.
#'
#' @param x See \code{\link[base]{Extract}}, mgnet object.
#' @param i See \code{\link[base]{Extract}}, samples indices.
#' @param j See \code{\link[base]{Extract}}, taxa indices
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
setMethod(f="[",
          signature="mgnet",
          definition=function(x,i,j){
            
            ifelse(length(x@data)!=0, data.new<-x@data[i,j,drop=F], data.new<-x@data)
            ifelse(length(x@meta_sample)!=0, meta_sample.new<-x@meta_sample[i, ,drop=F], meta_sample.new<-x@meta_sample)
            ifelse(length(x@taxa)!=0, taxa.new<-x@taxa[j, ,drop=F], taxa.new<-x@taxa)
            ifelse(length(x@meta_taxa)!=0, meta_taxa.new<-x@meta_taxa[j, ,drop=F], meta_taxa.new<-x@meta_taxa)
            ifelse(length(x@log_data)!=0, log_data.new<-x@log_data[i,j,drop=F], log_data.new<-x@log_data)
            
            
            if(length(x@netw)==0){
              return(mgnet(data=data.new,
                           meta_sample=meta_sample.new,
                           taxa=taxa.new,
                           meta_taxa=meta_taxa.new,
                           log_data=log_data.new
              ))
            } else if(length(x@netw)!=0 & missing(i)){
              return(selection_taxa(x,j))
            } else {
              warning("Sample subsetting is not defined for network the resulting object will not have the slots netw and comm")
              return(mgnet(data=data.new,
                           meta_sample=meta_sample.new,
                           taxa=taxa.new,
                           meta_taxa=meta_taxa.new,
                           log_data=log_data.new
              ))
            }
          })
################################################################################
################################################################################
# END EXTRACTOR MG
################################################################################
################################################################################




################################################################################
################################################################################
# MGNET ARRANGE
################################################################################
################################################################################
#' Re-arranging the samples or taxa positioning.
#' 
#' @description The function permits to order differently the samples or the taxa
#' inside the mgnet object.
#' 
#' @param object mgnet.
#' @param i sample indices that are numeric or character vectors or empty (missing). 
#' @param j taxa indices that are numeric or character vectors or empty (missing). 
#' 
#' @importFrom igraph graph_from_adjacency_matrix
#' @rdname arrange_mgnet-methods
#' @docType methods
#' @export
setGeneric("arrange_mgnet", function(object,i,j) standardGeneric("arrange_mgnet"))
#' @rdname arrange_mgnet-methods
setMethod("arrange_mgnet", "mgnet",
          function(object,i,j){
            
            # checks
            #------------------------------------------------------------------#
            if(!missing(i)){
              
              if(!(is.character(i)|is.numeric(i)|is.function(i))) stop("i must be numeric or character or a function")
              if(is.function(i)) i <- i(object)
              if(!is.null(dim(i))) stop("i must must be a vector")
              if(length(i)!=nsample(object)) stop("i must have the length equal to the sample number of object")
              if(any(duplicated(i))) stop("find at least a duplicated in i")
              if(is.character(i)){
                if(any(!(i%in%sample_name(object)))) stop("find at least an element in i not present in the sample_name of object")
              }
              if(is.numeric(i)){
                if(any(!(i%in%(1:nsample(object))))) stop("i elements must be in range 1 to sample number of object")
              }
            }
            
            if(!missing(j)){
              
              # checks
              if(!(is.character(j)|is.numeric(j)|is.function(j))) stop("i must be numeric or character")
              if(is.function(j)) j <- j(object)
              if(!is.null(dim(j))) stop("j must must be a vector")
              if(length(j)!=ntaxa(object)) stop("i must have the length equal to the taxa number of object")
              if(any(duplicated(j))) stop("find at least a duplicated in j")
              if(is.character(j)){
                if(any(!(j%in%taxaID(object)))) stop("find at least an element in j not present in the taxaID of object")
              }
              if(is.numeric(j)){
                if(any(!(j%in%(1:ntaxa(object))))) stop("j elements must be in range 1 to taxa number of object")
              }
            }
            #----------------------------------------------------------------#
            # end checks
            
            if(missing(i)) i <- 1:nsample(object)
            if(missing(j)) j <- 1:ntaxa(object)
            
            ifelse(length(object@data)!=0, data.new <- object@data[i,j],
                                           data.new <- object@data)
            ifelse(length(object@meta_sample)!=0, meta_sample.new <- object@meta_sample[i,],
                                                  meta_sample.new <- object@meta_sample)
            ifelse(length(object@taxa)!=0, taxa.new <- object@taxa[j,],
                                           taxa.new <- object@taxa)  
            ifelse(length(object@meta_taxa)!=0, meta_taxa.new <- object@meta_taxa[j,],
                                                meta_taxa.new <- object@meta_taxa)
            ifelse(length(object@log_data)!=0, log_data.new <- object@log_data[i,j],
                                               log_data.new <- object@log_data)
            
            if(length(object@netw)!=0){
              adjM <- adjacency_matrix(object)
              adjM <- adjM[j,j]
              netw.new <- graph_from_adjacency_matrix(adjM, mode="undirected",weighted=T)
            } else {
              netw.new <- object@netw
            }
            
            if(length(object@comm)!=0){
              comm.new <- object@comm
              comm.new$membership <- comm.new$membership[j]
            } else {
              comm.new <- object@comm
            }
            
            return(mgnet(data=data.new, meta_sample=meta_sample.new,
                         taxa=taxa.new, meta_taxa=meta_taxa.new,
                         log_data=log_data.new,
                         netw=netw.new, comm=comm.new))
              
          })
#' @rdname arrange_mgnet-methods
setMethod("arrange_mgnet","list",
          function(object,i,j){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="arrange_mgnet",signature="mgnet"),
                   i=i, j=j)}
)
################################################################################
################################################################################
# END ARRANGE MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# SHOW METHOD MGNET
################################################################################
################################################################################
#'@importFrom igraph vcount ecount edge_density membership sizes E
setMethod("show","mgnet",
          function(object){
            cat("******* Class mgnet , method Show ******* \n")
            
            if(!is.null(nsample(object))){
              cat(paste("Sample Number:",nsample(object),"\n"))
            }
            
            if(!is.null(nsample(object))){
              cat(paste("Taxa Number:",ntaxa(object),"\n"))
            }
            
            if(length(object@data)!=0){
              cat(paste("Zeros Percentage: ~",100*round(sum(object@data==0)/(nrow(object@data)*ncol(object@data)),4),"%\n",sep=""))
            }
            
            if(length(object@meta_sample)!=0){
              cat(paste("Sample Meta Info:",paste(colnames(object@meta_sample),collapse="," )),"\n")
            }
            
            if(length(object@taxa)!=0){
              cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")),"\n")
            }
            
            if(length(object@meta_taxa)!=0){
              cat(paste("Taxa Meta Info:",paste(colnames(object@meta_taxa),collapse="," )),"\n")
            }
            
            if(length(object@netw)!=0){
              cat(paste("Edge Number (Density): ",ecount(object@netw)," (~",round(igraph::edge_density(object@netw)*100,2),"%)","\n",sep=""))
              cat(paste("Positive Edge: ",sum(E(netw(object))$weight>0)," (~",100*round(sum(E(netw(object))$weight>0)/ecount(netw(object)),2),"%)\n",sep=""))
            }
            
            if(length(object@comm)!=0){
              
              if(length(object@comm)!=0){
                
                if("0" %in% names(sizes(object@comm))){
                  cat(paste("Signed Communities Number:",length(sizes(object@comm))-1,"\n"))
                  cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=","),"\n"))
                  cat(paste("Isolated Nodes:", sizes(object@comm)[[1]],"\n"))
                } else {
                  cat(paste("Signed Communities Number:",length(sizes(object@comm)),"\n"))
                  cat(paste("Communities Sizes:",paste(sizes(object@comm),collapse=","),"\n"))
                  cat("There aren't isolated nodes \n")
                }
              }
            }
            
            if(is.null(nsample(object)) & is.null(ntaxa(object))){
              cat("EMPTY\n")
            }
            
            cat("*********** End Show (mgnet) *********** \n")
          })
################################################################################
################################################################################
# END SHOW METHOD MGNET
################################################################################
################################################################################



################################################################################
################################################################################
# MGMELT
################################################################################
################################################################################
#' Melt mg object into form suitable for easy casting.
#' 
#' @description 
#' Summarize all elements present in an mg object in a single data frame. 
#' The functioning is similar to melt function of reshape2 package and it will 
#' be useful as a preprocess for ggplot2.
#' 
#' @usage mgmelt(object)
#' 
#' @param object mgnet.
#' 
#' @importFrom reshape2 melt
#' @rdname mgmelt-methods
#' @docType methods
#' @export
setGeneric("mgmelt", function(object) standardGeneric("mgmelt"))
#' @rdname mgmelt-methods
setMethod("mgmelt", "mgnet",
          function(object){
            
            if(length(object@data)==0){stop("\n data slot must be present")}
            
            mdf <- reshape2::melt(data=object@data)
            colnames(mdf) <- c("SampleID","TaxaID","Abundance")
            rownames(mdf) <- paste(mdf$SampleID,"-",mdf$TaxaID,sep="")
            
            if("sample_sum"%in%sample_info(object)){
              mdf$Relative <- reshape2::melt(relative(object))$value
            }
            
            if(length(object@log_data)!=0){
              mdf$log_data <- reshape2::melt(object@log_data)$value
            }
            
            if(length(object@taxa!=0)) mdf <- cbind(mdf,object@taxa[mdf$TaxaID,])
            
            if(length(object@meta_sample!=0)){
              if(ncol(object@meta_sample)==1){
                df.tmp <- as.data.frame(object@meta_sample[mdf$SampleID,])
                colnames(df.tmp) <- colnames(object@meta_sample)
                mdf <- cbind(mdf,df.tmp)
              } else {
                mdf <- cbind(mdf, object@meta_sample[mdf$SampleID,])
              }
            } 
            
            if(length(object@meta_taxa!=0)){
              if(ncol(object@meta_taxa)==1){
                df.tmp <- as.data.frame(object@meta_taxa[mdf$TaxaID,])
                colnames(df.tmp) <- colnames(object@meta_taxa)
                mdf <- cbind(mdf,df.tmp)
              } else {
                mdf <- cbind(mdf, object@meta_taxa[mdf$TaxaID,])
              }
            } 
            
            if(length(object@comm)!=0) mdf <- cbind(mdf, "commID"=commID(object)[mdf$TaxaID])
            
            if(any(duplicated(t(mdf)))){
              mdf <- mdf[,-which(duplicated(t(mdf)))]
            }
            
            return(mdf)
          })
#' @rdname mgmelt-methods
setMethod("mgmelt","list",
          function(object){
            is_mgnet_list(object)
            res <- lapply(object, selectMethod(f="mgmelt",signature="mgnet"))
            res <- lapply(names(res), function(x) res[[x]] <- cbind("ID"=x,res[[x]]))
            res <- do.call("rbind", res)
            return(res)
          })
################################################################################
################################################################################




################################################################################
################################################################################
# SIGNED LAYOUT
################################################################################
################################################################################
#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.
#' 
#' @param obj network belong to igraph or mgnet class.
#'
#' @importFrom igraph layout.fruchterman.reingold subgraph.edges E is.igraph is_weighted
#' @rdname layout_signed-methods
#' @docType methods
#' @export
setGeneric("layout_signed", function(obj) standardGeneric("layout_signed"))
#' @rdname layout_signed-methods
setMethod("layout_signed","igraph",function(obj){
  
  if(!is_weighted(obj)) stop("obj must be a weighted graph")
  graph.sub <- subgraph.edges(graph=obj,
                              eids=which(E(obj)$weight>0),
                              delete.vertices=FALSE)
  
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
})
#' @rdname layout_signed-methods
setMethod("layout_signed","mgnet",function(obj){
  
  graph.sub <- subgraph.edges(graph=netw(obj),
                              eids=which(E(netw(obj))$weight>0),
                              delete.vertices=FALSE)
  
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
})
################################################################################
################################################################################




################################################################################
################################################################################
# REMOVE SMALLER COMMUNITIES
################################################################################
################################################################################
#' Remove smaller communities
#' 
#' @description Allows you to remove communities based on the number of vertices.
#' 
#' @param obj mgnet.
#' @param size integer indicates the vertex number threshold
#' @param trim (Optional, default TRUE) logical. If true, the function removes 
#' all nodes not belonging to a community with size equal to size. If false the 
#' filtered vertices are set as isolated.
#'  
#' @importFrom igraph vcount induced.subgraph V
#' @rdname remove_smaller_comm-methods
#' @docType methods
#' @export
setGeneric("remove_smaller_comm", function(obj,size,trim=TRUE) standardGeneric("remove_smaller_comm"))
#' @rdname remove_smaller_comm-methods
setMethod("remove_smaller_comm", c("mgnet","numeric"), function(obj, size, trim=TRUE){
  
  #Checks Arguments
  if(!is.numeric(size) | round(size)!=size | size<0) stop("size must be integer greater than 0")
  if(!is.logical(trim)) stop("keep must be logical")
  if(length(obj)==0) stop("missing communities info")
  #End Checks
  
  graph <- obj@netw
  comm <- obj@comm
  
  keep.comm.names <- as.numeric(names(sizes(comm)[sizes(comm)>=size]))
  keep.comm.vids <-  which(comm$membership %in% keep.comm.names)
  
  if(!trim){
    data <- obj@data
    taxa <- obj@taxa
    log_data <- obj@log_data
    meta_taxa <- obj@meta_taxa
    graph.sub <- graph
    comm.sub <- comm
    comm.sub$membership[setdiff(1:comm$vcount, keep.comm.vids)] <- 0
    comm.sub$modularity <- NA
  } else {
    data <- obj@data[,keep.comm.vids]
    taxa <- obj@taxa[keep.comm.vids,]
    log_data <- obj@log_data[,keep.comm.vids]
    meta_taxa <- obj@meta_taxa[keep.comm.vids,]
    graph.sub <- induced.subgraph(graph, V(graph)[keep.comm.vids])
    comm.sub <- comm
    comm.sub$membership <- comm$membership[keep.comm.vids]
    comm.sub$vcount <- length(comm.sub$membership)
    comm.sub$modularity <- NA
  }
  
  
  return(new("mgnet",
             data=data, meta_sample=obj@meta_sample, taxa=taxa, 
             meta_taxa=meta_taxa, netw=graph.sub, comm=comm.sub,
             log_data=log_data))
})
#' @rdname remove_smaller_comm-methods
setMethod("remove_smaller_comm",c("list","numeric"),
          function(obj,size,trim=TRUE){
            is_mgnet_list(obj)
            lapply(obj, selectMethod(f="remove_smaller_comm",
                                     signature=c("mgnet","numeric")),
                   size=size,trim=trim)})
################################################################################
################################################################################
# END REMOVE SMALLER COMMUNITIES
################################################################################
################################################################################




################################################################################
################################################################################
# ISOLATE 
################################################################################
################################################################################
#' Isolate Vertices
#' 
#' @description Remove vertices from their communities.
#' 
#' @usage isolate(object, ..., condition="AND")
#' 
#' @param object mgnet.
#' @param ... arbitrary set of vectors of functions. If vectors they must contain
#' numeric or name or logical indices to be preserved. If functions they must return
#' indices as just explained.
#' @param condition (Optional default "AND") character indicates the logical operation between indices.
#' 
#' @importFrom igraph make_clusters
#' @rdname isolate-methods
#' @docType methods
#' @export
setGeneric("isolate", function(object,...,condition="AND") standardGeneric("isolate"))
#' @rdname isolate-methods
setMethod("isolate", "mgnet", function(object,...,condition="AND"){
  
  condition <- match.arg(condition, c("AND","OR"))
  if(length(object@netw)==0) stop("netw slot cannot be empty")
  
  
  IDX <- list(...)
  for(i in 1:length(IDX)){
    
    idx <- IDX[[i]]
    if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
    
    if(is.function(idx)){
      idx <- idx(object)
      if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
    }  
    
    
    if(is.numeric(idx)){
      #numeric
      if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>ntaxa(object)){
        stop("numeric indices must be integers in range 1 to the maximum nuber of taxa")
      }
      #transform to logical
      IDX[[i]] <- (1:nsample(object))%in%idx
      
    }else if(is.character(idx)){
      #character
      if( !any(idx%in%taxa_name(object)) ){
        stop("string indices must be a subset of taxa_name of object")
      }
      #transform to logical
      IDX[[i]] <- taxa_name(object)%in%idx
      
    }else if(is.logical(idx)){
      #logical
      if(length(idx)!=ntaxa(object)){
        stop("logical indices must have the length equal to taxa number")
      }
      IDX[[i]] <- idx
      
    } else {stop("the indices must character, numeric or logical vectors")}
  }
  
  
  IDX <- do.call(rbind,IDX)
  # Condition  
  if(condition=="AND"){
    IDX <- as.logical(apply(IDX,2,prod))
  } else {
    IDX <- apply(IDX,2,sum)>0
  }
  
  membership <- comm(object)$membership
  membership[IDX] <- 0
  
  comm.new <- make_clusters(netw(object),
                            membership=membership,
                            algorithm="modified",
                            modularity=1)
  comm.new$modularity <- NA
  
  return(mgnet(data=object@data,
               meta_sample=object@meta_sample,
               taxa=object@taxa,
               meta_taxa=object@meta_taxa,
               log_data=object@log_data,
               netw=object@netw,
               comm=comm.new))
  
})
#' @rdname isolate-methods
setMethod("isolate","list",
          function(object,...,condition="AND"){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="isolate",
                                        signature="mgnet"),
                   ...=..., condition=condition)}
)
################################################################################
################################################################################
# END ISOLATE
################################################################################
################################################################################




################################################################################
################################################################################
# DEGREE MGNET
################################################################################
################################################################################
#' Signed Vertices Degrees with Communities Info
#'
#' @description Calculates the vertices degree using the information on the 
#' sign of the edges and the communities. The function make possible to distinguish
#' between positive, negative or both sign of the edges and discern intra or extra
#' (also both cases) communities edges.
#' 
#' @param obj mgnet class.
#' @param sign (default "all") character indicates the sign of the edges. Possible values are
#' "positive","negative","all".
#' @param type (default "all") character with possible values "intra","extra","all".
#'
#' @importFrom stats setNames
#' @importFrom igraph degree subgraph.edges intersection crossing
#' @rdname degree_mgnet-methods
#' @docType methods
#' @export
setGeneric("degree_mgnet", function(obj,sign="all",type="all") standardGeneric("degree_mgnet"))
#' @rdname degree_mgnet-methods
setMethod("degree_mgnet","mgnet",function(obj,sign="all",type="all"){
  
  sign <- match.arg(sign,c("positive","negative","all"))
  type <- match.arg(type,c("intra","extra","all"))
  
  if(type!="all" & length(comm(obj))==0) stop("With type setted also comm slots must be no empty")
  
  n <- netw(obj)
  if(length(obj@comm)!=0) c <- comm(obj)

  
  if(sign=="positive"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight>0],delete.vertices=FALSE)
  } else if(sign=="negative"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight<0],delete.vertices=FALSE)
  } else {
    sub.sign <- n
  }
  
  # Modify isolated nodes as igraph need
  #------------------------------------#
  if(type!="all"){
    c$membership[c$membership==0] <- (length(c)+1):((length(c)+1)+length(c$membership[c$membership==0])-1)
  }
  
  if(type=="intra"){
    sub.type <- subgraph.edges(graph=n,eids=E(n)[!crossing(c,n)],delete.vertices=FALSE)
  } else if(type=="extra"){
    sub.type <- subgraph.edges(graph=n,eids=E(n)[crossing(c,n)],delete.vertices=FALSE)
  } else {
    sub.type <- n
  }
  
  g.int <- intersection(sub.sign,sub.type)
  return(degree(g.int))
})
#' @rdname degree_mgnet-methods
setMethod("degree_mgnet","list",
          function(obj,sign="all",type="all"){
            is_mgnet_list(obj)
            lapply(obj, selectMethod(f="degree_mgnet",
                                     signature="mgnet"),
                   sign=sign,type=type)})
################################################################################
################################################################################
# END DEGREE MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# STRENGTH MGNET
################################################################################
################################################################################
#' Signed Vertices Strengths with Communities Info
#'
#' @description Calculates the vertices strengths using the information on the 
#' sign of the edges and the communities. The function make possible to distinguish
#' between positive, negative or both sign of the edges and discern intra or extra
#' (also both cases) communities edges.
#' 
#' @param obj mgnet class.
#' @param sign (default "all") character indicates the sign of the edges. Possible values are
#' "positive","negative","all".
#' @param type (default "all") character with possible values "intra","extra","all".
#'
#' @importFrom igraph degree subgraph.edges intersection
#' @rdname strength_mgnet-methods
#' @docType methods
#' @export
setGeneric("strength_mgnet", function(obj,sign="all",type="all") standardGeneric("strength_mgnet"))
#' @rdname strength_mgnet-methods
setMethod("strength_mgnet","mgnet",function(obj,sign="all",type="all"){
  
  sign <- match.arg(sign,c("positive","negative","all"))
  type <- match.arg(type,c("intra","extra","all"))
  
  if(type!="all" & length(comm(obj))==0) stop("With type setted also comm slots must be no empty")
  
  n <- netw(obj)
  if(length(obj@comm)!=0) c <- comm(obj)
  
  # Modify isolated nodes as igraph want
  #------------------------------------#
  c$membership[c$membership==0] <- (length(c)+1):((length(c)+1)+length(c$membership[c$membership==0])-1)
  
  if(sign=="positive"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight>0],delete.vertices=FALSE)
  } else if(sign=="negative"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight<0],delete.vertices=FALSE)
  } else {
    sub.sign <- n
  }
  
  if(type=="intra"){
    sub.type <- subgraph.edges(graph=n,eids=E(n)[!crossing(c,n)],delete.vertices=FALSE)
  } else if(type=="extra"){
    sub.type <- subgraph.edges(graph=n,eids=E(n)[crossing(c,n)],delete.vertices=FALSE)
  } else {
    sub.type <- n
  }
  
  g.int <- intersection(sub.sign,sub.type)
  adj.int <- as_adjacency_matrix(g.int,attr="weight_1",sparse=FALSE)
  strength <- setNames(rowSums(adj.int),taxaID(obj))
  return(strength)
})
#' @rdname strength_mgnet-methods
setMethod("strength_mgnet","list",
          function(obj,sign="all",type="all"){
            is_mgnet_list(obj)
            lapply(obj, selectMethod(f="strength_mgnet",
                                     signature="mgnet"),
                   sign=sign,type=type)})
################################################################################
################################################################################
# END STRENGTH MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# INTRA BETWEENESS MGNET
################################################################################
################################################################################
#' Inta Communities Betweenness Centrality 
#'
#' @description Calculates the betweenness centralities of each vertices inside
#' its communities. For the scope the function take into account only positive 
#' links.
#' 
#' @param obj mgnet class.
#' @param normalize allows to normalize the betweenness centralities scores.
#' "nVertices" if it will be divided respect the vertex number, "L1" if the sum
#' constrain is equal to 1 and "none" to leave invariate.
#'
#' @importFrom igraph betweenness delete_vertices membership delete_edges
#' @rdname betweenness_mgnet-methods
#' @docType methods
#' @export
setGeneric("betweenness_mgnet", function(obj, normalize="none") standardGeneric("betweenness_mgnet"))
#' @rdname betweenness_mgnet-methods
setMethod("betweenness_mgnet","mgnet",function(obj, normalize="none"){
  
  normalize <- match.arg(normalize,c("none","nVertices","L1"))
  if(length(comm(obj))==0) stop("comm slot cannot be empty")
  
  btw <- setNames(rep(0,ntaxa(obj)), taxaID(obj))
  
  for(c in 1:max(membership(comm(obj))) ){
    
    sub <- delete_vertices(netw(obj), membership(comm(obj))!=c )
    sub <- delete_edges(sub, which(E(sub)$weight<0))
    sub <- delete_vertices(sub, degree(sub)==0)
    
    if(normalize=="none"){
      btw.tmp <- betweenness(sub, directed=FALSE, normalized=FALSE)
    } else if (normalize=="nVertices") {
      btw.tmp <- betweenness(sub, directed=FALSE, normalized=TRUE)
    } else {
      btw.tmp <- betweenness(sub, directed=FALSE, normalized=FALSE)
      btw.tmp <- btw.tmp/sum(btw.tmp)
    }
    
    btw[names(btw.tmp)] <- btw.tmp
  }
  
  return(btw)
})
#' @rdname betweenness_mgnet-methods
setMethod("betweenness_mgnet","list",
          function(obj,normalize="none"){
            is_mgnet_list(obj)
            lapply(obj, selectMethod(f="betweenness_mgnet",
                                     signature="mgnet"),
                   normalize=normalize)})
################################################################################
################################################################################
# END INTRA BETWEENNESS MGNET
################################################################################
################################################################################




################################################################################
################################################################################
# PLOT MGNET
################################################################################
################################################################################
#' plot.mgnet
#' 
#' @description
#' Plot network belonging to a mgnet object. The function focus on characterize
#' the signed properties of the edges, color them and placing the vertices considering
#' only the positive ones. Furthermore it permits to tune some graphical aspects
#' of the nodes and of the edges. At vertex.size, named \eqn{v_0}, of the graph is applied
#' the following function:
#' \deqn{v=\text{\footnotesize{sumConst}}+\text{\footnotesize{multConst}}\left(\text{\footnotesize{max}}(v_0)\left(\frac{v_0}{\text{\footnotesize{max}}(v_0)}\right)^{\text{expFactor}}\right)}
#' 
#' @param x mgnet object
#' @param layout vertices layout in igraph style. If missing it is applied the layout_signed.
#' @param vertexSize numeric or character indicated the size of the vertices. If 
#' is character the variable can assumes the two values 'log-mean' or 'var-mean' and
#' where the vertices sizes are taken from the colMeans or from the columns variance of
#' the log_data slot. Otherwise if vertexSize is numeric the vertices scales as received.
#' @param expFactor numeric value and It tune the difference in size between the smaller and the higher sizes of vertices.
#' @param multConst numeric positive applied to all the vertices size as multiplicative factor.
#' @param sumConst numeric positive summed to all vertices size.
#' @param thickness logical that it scale the width of the edges using the its weights in absolute values.
#' @param alphaFactor numeric in range \[0,1\] and it modify the alpha trasparency of edges color.
#' @param widthFactor numeric positive value that scale the width of all edges.
#' @param posCol rgb color for positive edges.
#' @param negCol rgb color for negative edges.
#' @param maxSize numeric positive integer with default NULL value. If it is set
#' re-scale all the vertices size in order to obtain the larger vertices with the
#' size equal to maxSize.
#' @param maxWidth numeric positive integer with default NULL value. If it is set
#' re-scale all the edges width in order to obtain the larger edge with the
#' width equal to maxWidth.
#' @param ... additional arguments to be passed to igraph plot of the network
#' 
#' @importFrom igraph E<- V<-
#' @importFrom stats var
#' @importFrom grDevices adjustcolor
#' @rdname plot.mgnet
#' @export
plot.mgnet <- function(x,layout,
                       vertexSize="mean-log",
                       expFactor=1, multConst=1, sumConst=0,
                       thickness=TRUE,
                       alphaFactor=.5, widthFactor=1,
                       posCol=rgb(0,0,1), negCol=rgb(1,0,0),
                       maxSize=NULL, maxWidth=NULL,
                       ...) {
  if(length(netw(x))==0) stop("missing network in mgnet")
  
  if(is.character(vertexSize)){
    
    if(length(x@log_data)==0) stop("if vertexSize is equal to mean-log or var-log the slot log_data cannot be empty")
    vertexSize <- match.arg(vertexSize, c("mean-log","var-log"))
    
    if(vertexSize=="mean-log"){
      minV <- ifelse(any(x@log_data<0), 
                     abs(min(x@log_data)),
                     0)
      vertexSize <- colMeans(x@log_data + minV)
    } else {
      vertexSize <- apply(x@log_data,2,stats::var)
    }
    
  } else if (is.numeric(vertexSize)){
    if(length(vertexSize)==1){
      if(vertexSize<0) stop("vertexSize cannot be negative")
      vertexSize<-rep(vertexSize, ntaxa(x))
    } else if(length(vertexSize)==ntaxa(x)){
      if(any(vertexSize<0)) stop("vertexSize cannot contains negative values")
    }
  } else {
    stop("vertexSize can be 'mean-log', 'var-log' or a number or a numeric vector with all elelemnts greater than zero")
  }
  
  if(!is.numeric(expFactor)) stop("expFactor must be numeric")
  if(!is.numeric(multConst) | multConst<=0) stop("multFact must be a number greater than 0")
  if(!is.numeric(sumConst)) stop("sumConst must be a number")
  if(!is.logical(thickness)) stop("thickness must be logical")
  if(!is.numeric(alphaFactor) | alphaFactor<0 | alphaFactor>1) stop("alphaFactor must be a number in range [0,1]")
  if(!is.numeric(widthFactor) | widthFactor<=0) stop("widthFactor must be a positive number")
  
  if(missing(layout)){
    layout <- layout_signed(netw(x))
  }
  
  g <- netw(x)
  
  # Edges
  if(thickness){
    E(g)$width <- abs(E(g)$weight)
    E(g)$width <- widthFactor*E(g)$width
    if(!is.null(maxWidth)){
      E(g)$width <- (E(g)$width/max(E(g)$width))*maxWidth
    }
  } else {
    E(g)$width <- 1
  }
  E(g)$color <- ifelse(E(g)$weight>0,posCol,negCol)
  E(g)$color <- adjustcolor(E(g)$color,alpha.f=alphaFactor)
  
  # Nodes
  V(g)$size <- max(vertexSize)*(vertexSize/max(vertexSize))^(expFactor)
  V(g)$size <- sumConst + multConst* V(g)$size
  if(!is.null(maxSize)){
    V(g)$size <- (V(g)$size/max(V(g)$size))*maxSize
  }
  
  plot(g,layout=layout,...)
}
################################################################################
################################################################################
# END PLOT MGNET
################################################################################
################################################################################



