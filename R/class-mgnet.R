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
      if(any(sapply(colnames(object@data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("taxa/col name cannot begin with number")
      if(any(sapply(rownames(object@data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("sample/row name cannot begin with number")
    }
    
    #CHECK META_SAMPLE
    #-------------------------------------#
    if( length(object@meta_sample)!=0 ){
      if(is.null(rownames(object@meta_sample))) return("\n meta_sample data.frame must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@meta_sample))) return("\n meta_sample data.frame have the cols names where the experimental variables were saved")
      if(any(duplicated(rownames(object@meta_sample)))) return("\n find in meta_sample matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@meta_sample)))) return("\n find in meta_sample matrix at least a duplicated col name / experimental variable.")
      if(any(sapply(rownames(object@data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("sample/row name cannot begin with number")
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
      if(any(sapply(colnames(object@data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("taxa/col name cannot begin with number")
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
      if(any(sapply(rownames(object@meta_taxa),function(x)grepl(substr(x,1,1),'[0-9]')))) return("taxa/row name in meta_taxa cannot begin with number")
    }
    
    #CHECK LOG_DATA
    #-------------------------------------#
    if( length(object@log_data)!=0 ){
      if(!is.numeric(object@log_data)) return("\n log_data matrix must be numeric")
      if(is.null(rownames(object@log_data))) return("\n log_data matrix must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@log_data))) return("\n log_data matrix must have the cols names where the taxa IDs were saved")
      if(any(duplicated(rownames(object@log_data)))) return("\n find in log_data matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@log_data)))) return("\n find in log_data matrix at least a duplicated col name / taxa ID.")
      if(any(sapply(colnames(object@log_data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("taxa/col name cannot begin with number")
      if(any(sapply(rownames(object@log_data),function(x)grepl(substr(x,1,1),'[0-9]')))) return("sample/row name cannot begin with number")
    }
    
    # CHECK NETW
    #-------------------------------------#
    if(length(object@netw)!=0){
      if(is_directed(object@netw)) return("\n netw slot must be an undirected igraph object")
      if(!is_weighted(object@netw)) return("\n netw slot must be an weighted igraph object")
      if(is.null(V(object@netw)$name)) return("\n netw vertices name ( V(netw(obj))$name ) cannot be empty.")
      if(any(sapply(V(object@netw)$name,function(x)grepl(substr(x,1,1),'[0-9]')))) return("taxa name in netw cannot begin with number")
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
  check <- sapply(L, function(x) isa(x,"mgnet"))
  
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
  
  if(length(object@meta_sample)!=0){
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
#' @param adj adjacency matrix of netw.
#' @param netw un-directed weighted igraph network.
#' @param comm communities object associated to netw.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
mgnet <- function(data=matrix(nrow=0,ncol=0),
                  meta_sample=data.frame(),
                  taxa=matrix(nrow=0,ncol=0),
                  meta_taxa=data.frame(),
                  adj=matrix(nrow=0,ncol=0),
                  log_data=matrix(nrow=0,ncol=0),
                  netw=make_empty_graph(n=0, directed=FALSE),
                  comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))
                  ){
  
  if(length(netw)!=0 & length(adj)!=0){
    stop("the 'adj' and 'netw' arguments cannot be specified together")
  }
  
  if(length(adj)!=0){
    if(!is.numeric(adj) | !is.matrix(adj) | !isSymmetric(adj)) stop("adj must be simmetric matrix")
  }
  
  if(length(netw)==0 & length(adj)!=0){
    netw <- graph_from_adjacency_matrix(adj,'undirected',weighted=TRUE)
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
#'
#' @export
#' @docType methods
#' @rdname save_log_data-methods
setGeneric("save_log_data", function(object, method="base") standardGeneric("save_log_data"))
#' @rdname save_sample_sum-methods
setMethod("save_log_data", "mgnet", function(object, method="base"){
  
  method <- match.arg(method, c("base","inter-quantile"))
  
  ifelse(any(data(object)==0),
         x <- zero_dealing(data(object)),
         x <- data(object))
  
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
setMethod("log_data<-", c("mgnet", "data.frame"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=object@taxa, log_data=value,
      log_data=object@log_data,
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
#' @param value \code{\link{igraph}}
#'
#' @export
#' @docType methods
#' @rdname assign-netw
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))
#' @rdname assign-netw
setMethod("netw<-", c("mgnet", "igraph"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=object@taxa, meta_taxa=object@meta_taxa,
      log_data=object@log_data,
      netw=value, comm=object@comm)
})
#####################################
# COMM<-
#####################################
#' Assign a new comm to \code{object}
#'
#' @usage comm(object) <- value
#'
#' @param object mgnet.
#' @param value \code{\link{communities}}
#'
#' @export
#' @docType methods
#' @rdname assign-comm
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))
#' @rdname assign-comm
setMethod("comm<-", c("mgnet", "communities"), function(object, value){
  new("mgnet",data=object@data, meta_sample=object@meta_sample, 
      taxa=object@taxa, meta_taxa=object@meta_taxa,
      log_data=object@log_data,
      netw=object@netw, comm=value)
})
################################################################################
################################################################################
# END SETTERS MGNET
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
            cat(paste("Sample Number:",max(nrow(object@data),nrow(object@meta_sample)),"\n"))
            cat(paste("Taxa Number:",max(ncol(object@data),nrow(object@taxa),
                                         nrow(object@meta_taxa),vcount(object@netw)),"\n"))
            cat(paste("Zeros Percentage: ~",
                      100*round(sum(object@data==0)/(nrow(object@data)*ncol(object@data)),4),
                      "%\n",sep=""))
            cat(paste("Sample Meta Info:",paste(colnames(object@meta_sample),collapse="," )),"\n")
            cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")),"\n")
            cat(paste("Taxa Meta Info:",paste(colnames(object@meta_taxa),collapse="," )),"\n")
            cat(paste("Edge Number (Density): ",ecount(object@netw)," (~",
                      round(edge_density(object@netw)*100,2),"%)","\n",sep=""))
            cat(paste("Positive Edge: ",sum(E(netw(object))$weight>0)," (~",
                      100*round(sum(E(netw(object))$weight>0)/ecount(netw(object)),2),
                      "%)\n",sep=""))
            
            if(length(object@comm)!=0){
              
              if("0" %in% names(sizes(object@comm))){
                cat(paste("Signed Communities Number:",length(sizes(object@comm))-1,"\n"))
                cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=","),"\n"))
                cat(paste("Isolated Nodes:", sizes(object@comm)[[1]],"\n"))
              } else {
                cat(paste("Signed Communities Number:",length(sizes(object@comm)),"\n"))
                cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=","),"\n"))
                cat("There aren't isolated nodes \n")
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
  if(length(object@data!=0)) return(ncol(object@data))
  else if(length(object@taxa!=0)) return(nrow(object@taxa))
  else if(length(object@meta_taxa!=0)) return(nrow(object@meta_taxa))
  else if(length(object@log_data!=0)) return(ncol(object@log_data))
  else if(length(object@netw!=0)) return(vcount(object@netw))
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
setMethod("sample_sum", "mgnet",function(object)return(object@meta_sample$sample_sum))
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
#' @param object mgnet or list.
#' @param rank (Optional)
#' See \code{\link{save_sample_sum}}.
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
# CLR ABUNDANCE 
#####################################
#' Get clr abundances at taxonomy rank choosen.
#' 
#' @description 
#' Retrieves the clr abundances of data the taxonomic rank choosen.
#' In mgnet object must be save log_data slot.
#' 
#' @usage relative(object,rank)
#' 
#' @param object mgnet or list.
#' @param rank (Optional)
#' 
#' @importFrom stats aggregate
#' 
#' @rdname clr_abundance-methods
#' @docType methods
#' @export
setGeneric("clr_abundance", function(object,rank) standardGeneric("clr_abundance"))
#' @rdname clr_abundance-methods
setMethod("clr_abundance", c("mgnet","missing"),function(object){
  
  if(length(object@log_data)==0) stop("log_data cannot be empty")
  return(object@log_data)
  
})
#' @rdname clr_abundance-methods
setMethod("clr_abundance", c("mgnet","character"),function(object,rank){
  if(length(object@log_data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
  if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))

  different.taxa <- unique(object@taxa[,rank])
  data.aggregate <- t(aggregate(t(object@log_data), by=list(object@taxa[,rank]), 
                                FUN="sum", drop=FALSE))
  colnames(data.aggregate) <- data.aggregate[1,]
  data.aggregate <- data.aggregate[-1,]
  data.aggregate <- data.aggregate[,different.taxa]
  class(data.aggregate) <- "numeric"
  
  return(data.aggregate)
})
#' @rdname clr_abundance-methods
setMethod("clr_abundance",c("list","missing"),
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="clr_abundance",signature=c("mgnet","missing")))})
#' @rdname clr_abundance-methods
setMethod("clr_abundance",c("list","character"),
          function(object,rank){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="clr_abundance",signature=c("mgnet","character")),
                   rank=rank)}
)
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
#' @param object mgnet-class.
#'
#' @importFrom igraph membership
#' @rdname commID-methods
#' @docType methods
#' @export
setGeneric("commID", function(object) standardGeneric("commID"))
#' @rdname commID-methods
setMethod("commID", "mgnet", function(object){
  if(length(object@comm)!=0){
    return(as.character(membership(object@comm)))
  } else {
    stop("the comm slot is not present")
  }
})
#' @rdname commID-methods
setMethod("commID","list",
          function(object){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="commID",signature="mgnet"))})
#####################################
# TOP TAXA
#####################################
#' Get top taxa.
#' 
#' @description 
#' Sort the top taxa at taxonomic level and function choosen
#' 
#' @param object mgnet or list.
#' @param which the comparison can be made on absolute values, on relative abundances
#' or on clr transformed data (possible values are "abs","rel","clr").
#' @param rank taxonomic level choosen.
#' @param fun function choosen for the comparison. It must receve as input a vector
#' and return a number (mean, median, max, etc...).
#' 
#' 
#' @rdname top_taxa-methods
#' @docType methods
#' @export
setGeneric("top_taxa", function(object,which,rank,fun) standardGeneric("top_taxa"))
#' @rdname top_taxa-methods
setMethod("top_taxa", c("mgnet","character","character","function"),
          function(object,which,rank,fun){
            
            which <- match.arg(which,c("abs","rel","clr"))
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
            if(!is.function(fun)) stop("fun must be a function")
            if(!is.numeric(fun(c(1,2)))) stop("fun must return a number")
            if(length(fun(c(1,2)))!=1) stop("fun must return an element of length 1")
            
            if(which=="abs"){
              x <- abundance(object,rank)
            } else if(which=="rel"){
              if(!("sample_sum"%in%sample_info(object))) stop("for relative must be present sample_sum in meta_sample (See save_sample_sum).")
              x <- relative(object,rank)
            } else {
              if(length(log_data(object))==0) stop("log_data slot cannot be empty")
              x <- clr_abundance(object,rank)
            }
            
            top_taxa <- apply(x,2,FUN=fun)
            names(top_taxa) <- colnames(x)
            top_taxa <- sort(top_taxa,decreasing=TRUE)
            return(top_taxa)
  })
#' @rdname top_taxa-methods
setMethod("top_taxa","list",
          function(object,which,rank,fun){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="top_taxa",
                                        signature=c("mgnet","character","character","function")),
                   which=which,rank=rank,fun=fun)})
################################################################################
################################################################################
# END BASE METHODS
################################################################################
################################################################################




################################################################################
################################################################################
# AGGREGATE TAXA
################################################################################
################################################################################
#' Organize data in higher taxonomic level.
#' 
#' @description 
#' Reorganize an \code{\link{mgnet-class}} object in an higher taxonomy rank.
#' The function sums the taxa with the same classification and return a new mg
#' object with appropriate data and taxa slots.
#' 
#' @usage aggregate_taxa(object, rank)
#' 
#' @param object mgnet.
#' @param rank taxonomic level choosen.
#' 
#' @importFrom stats aggregate
#' 
#' @rdname aggregate_taxa-methods
#' @docType methods
#' @export
setGeneric("aggregate_taxa", function(object, rank) standardGeneric("aggregate_taxa"))
#' @rdname aggregate_taxa-methods
setMethod("aggregate_taxa", c("mgnet","character"),
          function(object, rank){
            
            if(length(object@data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
            
            different.taxa <- aggregate(object@data[1,],list(object@taxa[,rank]),sum)$Group.1
            data.aggregate <- t(apply(object@data,1,
                                    function(x) aggregate(x,list(object@taxa[,rank]),sum, drop=F)$x,
                                    simplify=T))
            colnames(data.aggregate) <- different.taxa
            
            taxa.aggregate <- object@taxa[,1:which(ranks(object)==rank)]
            taxa.aggregate <- taxa.aggregate[!duplicated(taxa.aggregate), ]
            rownames(taxa.aggregate) <- taxa.aggregate[,rank]
            data.aggregate <- data.aggregate[,rownames(taxa.aggregate)]
            
            return(new("mgnet",data=data.aggregate, meta_sample=object@meta_sample, taxa=taxa.aggregate))
          })
#' @rdname aggregate_taxa-methods
setMethod("aggregate_taxa",c("list","character"),
          function(object,rank){
            is_mgnet_list(object)
            lapply(object, selectMethod(f="aggregate_taxa",signature=c("mgnet","character")),
                   rank=rank)}
)
################################################################################
################################################################################
# END AGGREGATE TAXA
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
            
            if(length(log_data)!=0){
              mdf$CLR <- reshape2::melt(object@log_data)$value
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
            
            
            if(length(object@meta_taxa!=0)) mdf <- cbind(mdf,object@meta_taxa[mdf$TaxaID,])
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
# END MGMELT
################################################################################
################################################################################




################################################################################
################################################################################
# SAMPLE SELECTION
################################################################################
################################################################################
#' Selection a subset of samples from mgnet object.
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
                         log_data=log_data))
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
                           meta_taxa=meta_taxa.new
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
# ARRANGE VERTICES
################################################################################
################################################################################
#' Arrange the vertex number
#'
#'@description Edit the vertices/taxa of the \code{\link{mgnet}} object 
#'(for example to compare two different networks)
#'
#'@param obj \code{\link{mgnet}} object
#'@param new_taxa data.frame with taxonomic information about the vertices wanted in new graph.
#'  The dataset must have the same form as taxa in MetaGenomic/MetaGenomicGraph object.
#'  Naturally in new_taxa all the old vertices must be present.
#'@param new_meta_taxa (Optional) data.frame with additional info on taxa.
#'  
#'@importFrom igraph make_clusters
#'@rdname arrange_vertices-methods
#'@docType methods
#'@export
setGeneric(name="arrange_vertices",
           def=function(obj, new_taxa, new_meta_taxa=data.frame()) standardGeneric("arrange_vertices"))
#' @rdname arrange_vertices-methods
setMethod("arrange_vertices",c("mgnet","matrix"),
          function(obj, new_taxa, new_meta_taxa=data.frame()){
            
            if(is.null(rownames(new_taxa)) | is.null(colnames(new_taxa))) stop("new_taxa must have rownames (taxaID) and colnames (info) assigned")
            if(length(obj@taxa)==0 | length(obj@netw)==0) stop("taxa and netw cannot be empty")
            if(any(!(taxaID(obj)%in%rownames(new_taxa)))) stop("find at least a taxa not present in new_taxa")
            if(any(ranks(obj)!=colnames(new_taxa))) stop('new_taxa must have the same ranks as obj')
            if(length(new_meta_taxa)!=0){
              if(!is.data.frame(new_meta_taxa)) stop("new_meta_taxa must be a data.frame")
              if(is.null(rownames(new_meta_taxa)) | is.null(colnames(new_meta_taxa))) stop("new_meta_taxa must have rownames (taxaID) and colnames (info) assigned")
              if(any(rownames(new_taxa)!=rownames(new_meta_taxa))) stop("new_taxa and new_meta_taxa muse have the same rownames")
            }
            if(length(obj@meta_taxa)!=0 & length(new_meta_taxa)==0) stop("miss new_meta_taxa and on contrary meta_taxa of obj it is not empty")
            
            sample_name <- sample_name(obj)
            ntaxa <- nrow(new_taxa)
            taxa_name <- rownames(new_taxa)
            nsample <- nsample(obj)
            
            if(length(data)!=0){
              data <- matrix(0,nrow=nsample,ncol=ntaxa,
                             dimnames=list(sample_name,taxa_name))
              data[,taxaID(obj)] <- obj@data
            } else {
              data <- obj@data
            }

            if(length(log_data)!=0){
              log_data <- matrix(NA_real_,nrow=nsample,ncol=ntaxa,
                                 dimnames=list(sample_name,taxa_name))
              log_data[,taxaID(obj)] <- obj@log_data
            } else {
              log_data <- obj@log_data
            }

            adj <- as_adjacency_matrix(obj@netw,sparse=F,attr="weight")
            adj.new <- matrix(0,nrow=ntaxa,ncol=ntaxa,
                              dimnames=list(taxa_name,taxa_name))
            adj.new[taxaID(obj),taxaID(obj)] <- adj
            graph.new <- graph_from_adjacency_matrix(adj.new,mode="undirected",
                                                     weighted=T)
            
            if(length(obj@comm)!=0){
              
              mem <- membership(obj@comm)
              names(mem) <- taxaID(obj)
              
              mem.new <- rep(0, ntaxa)
              names(mem.new) <- taxa_name
              mem.new[names(mem)] <- mem
              
              comm.new <- make_clusters(graph=graph.new,
                                        membership=mem.new,
                                        algorithm="signed weights louvain",
                                        modularity=FALSE)
              comm.new$modularity <- NA
            } else {
              comm.new <- obj@comm
            }
            
            return(mgnet(data=data,meta_sample=obj@meta_sample,
                         taxa=new_taxa, meta_taxa=new_meta_taxa,
                         netw=graph.new,comm=comm.new,
                         log_data=log_data))
          })
#' @rdname arrange_vertices-methods
setMethod("arrange_vertices",c("list","matrix"),
          function(obj,new_taxa){
            is_mgnet_list(obj)
            lapply(obj, selectMethod(f="arrange_vertices",
                                     signature=c("mgnet","matrix")),
                   new_taxa=new_taxa)})
################################################################################
################################################################################
# END ARRANGE VERTICES
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
setMethod("remove_smaller_comm",c("list","numeric","logical"),
          function(obj,size,trim){
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
#' @importFrom igraph layout.fruchterman.reingold subgraph.edges E is.igraph
#' @rdname layout_signed-methods
#' @docType methods
#' @export
setGeneric("layout_signed", function(obj) standardGeneric("layout_signed"))
#' @rdname layout_signed-methods
setMethod("layout_signed","igraph",function(obj){
  
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
# END SIGNED LAYOUT
################################################################################
################################################################################




################################################################################
################################################################################
# PLOT MGNET
################################################################################
################################################################################
#'@importFrom igraph E<-
#'@export
plot.mgnet <- function(x,layout,...) {
  if(length(netw(x))==0) stop("missing network in mgnet")
  
  if(missing(layout)){
    layout <- layout_signed(netw(x))
  }
  
  E(netw(x))$color <- ifelse(E(netw(x))$weight>0, rgb(0,0,1,.5), rgb(1,0,0,.5))
  plot(netw(x),layout=layout,...)
}
################################################################################
################################################################################
# END PLOT MGNET
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