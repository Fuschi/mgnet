setOldClass("igraph")
setOldClass("communities")
################################################################################
################################################################################
# CLASS MGNET
################################################################################
################################################################################
#' S4 class to manage metagenomic networks
#'
#' @slot data table of data
#' @slot meta sample features
#' @slot taxa taxonomic information
#' @slot netw undirected, weighted and signed igraph network
#' @slot comm graph communities of netw
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy V vcount
#' @import methods
#' @name mgnet-class
#' @rdname mgnet-class
#' @exportClass mgnet
setClass(
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
# END CLASS MGNET
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
#' @param meta data.frame with experimental variables.
#' @param taxa character matrix with taxonomic classification.  
#' @param netw un-directed weighted igraph network.
#' @param comm communities object associated to netw.
#' @param adj adjacency matrix of netw.
#' @param mg object belong to mg class.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
mgnet <- function(data=matrix(nrow=0,ncol=0),
                  meta=data.frame(),
                  taxa=matrix(nrow=0,ncol=0),
                  netw=make_empty_graph(n=0, directed=FALSE),
                  comm=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE)),
                  adj=matrix(nrow=0,ncol=0),
                  mg=new("mg")){
  
  if((length(data)!=0 | length(meta)!=0 | length(taxa)!=0) & !empty(mg)){
    stop("if using an object of class mg then data, meta and taxa must be left empty")
  }
  
  if(length(netw)!=0 & length(adj)!=0){
    stop("the 'adj' and 'netw' arguments cannot be specified together")
  }
  
  if(length(adj)!=0){
    if(!is.numeric(adj) | !is.matrix(adj) | !isSymmetric(adj)) stop("adj must be simmetric matrix")
  }
  
  if(length(netw)==0 & length(adj)!=0){
    netw <- graph_from_adjacency_matrix(adj,'undirected',weighted=TRUE)
  }
  
  # if(length(netw)!=0 && length(comm)==0){
  #   comm <- cluster_signed(netw)
  # }
  
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
#' @rdname netw
#' @aliases netw,list
setMethod("netw","list",
          function(object){
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
#' @param object (Required) \code{\link{mgnet-class}}.
#'
#' @rdname comm
#' @docType methods
#' @export
setGeneric("comm", function(object) standardGeneric("comm"))
#' @rdname comm
#' @aliases comm,mgnet
setMethod("comm", "mgnet", function(object){return(object@comm)})
#' @rdname comm
#' @aliases comm,list
setMethod("comm","list",
          function(object){
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
          signature="mgnet",
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
setMethod("show","mgnet",
          function(object){
            cat("******* Class mgnet , method Show ******* \n")
            cat(paste("Sample Number:",max(nrow(object@data),nrow(object@meta)),"\n"))
            cat(paste("Taxa Number:",max(ncol(object@data),nrow(object@taxa)),"\n"))
            cat(paste("Zeros Percentage: ~",
                      100*round(sum(object@data==0)/(nrow(object@data)*ncol(object@data)),4),
                      "%\n",sep=""))
            cat(paste("Sample Meta Data:",paste(colnames(object@meta),collapse="," )),"\n")
            cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")),"\n")
            cat(paste("Edge Number (Density): ",ecount(object@netw)," (~",
                      round(edge_density(object@netw)*100,2),"%)","\n",sep=""))
            cat(paste("Positive Edge: ",sum(E(netw(object))$weight>0)," (~",
                      100*round(sum(E(netw(object))$weight>0)/ecount(netw(object)),2),
                      "%)\n",sep=""))
            
            if(length(object@comm)!=0){
              cat(paste("Signed Communities Number:",max(membership(object@comm)),"\n"))
              
              if("0" %in% names(sizes(object@comm))){
                cat(paste("Communities Sizes:",paste(sizes(object@comm)[-1],collapse=","),"\n"))
                cat(paste("Isolated Nodes:", sizes(object@comm)[[1]],"\n"))
              } else {
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
# NTAXA
#####################################
#' @importFrom igraph vcount
#' @rdname ntaxa
#' @aliases ntaxa,mgnet
setMethod("ntaxa", "mgnet", function(object){
  if(length(object@data!=0)) return(ncol(object@data))
  else if(length(object@taxa!=0)) return(nrow(object@taxa))
  else if(length(object@netw!=0)) return(vcount(object@netw))
  else return(NULL)
})
#' @rdname ntaxa
#' @aliases ntaxa,list
setMethod("ntaxa","list",
          function(object){
            lapply(object, selectMethod(f="ntaxa",signature="mgnet"))})
#####################################
# TAXAID
#####################################
#' @importFrom igraph V
#' @rdname taxaID
#' @aliases taxaID,mgnet
setMethod("taxaID", "mgnet", function(object){
  if(length(object@data!=0)) return(colnames(object@data))
  else if(length(object@taxa!=0)) return(rownames(object@taxa))
  else if(length(object@netw!=0)) return(V(object@netw)$name)
  else return(0)
})
#' @rdname taxaID
#' @aliases taxaID,list
setMethod("taxaID","list",
          function(object){
            lapply(object, selectMethod(f="taxaID",signature="mgnet"))})
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
#' @param object (Required) \code{\link{mgnet-class}}.
#'
#' @rdname commID
#' @docType methods
#' @export
setGeneric("commID", function(object) standardGeneric("commID"))
#' @importFrom igraph membership
#' @rdname commID
#' @aliases commID,mgnet
setMethod("commID", "mgnet", function(object){
  if(length(object@comm)!=0){
    return(as.character(membership(object@comm)))
  } else {
    stop("the comm slot is not present")
  }
})
#' @rdname commID
#' @aliases commID,list
setMethod("commID","list",
          function(object){
            lapply(object, selectMethod(f="commID",signature="mgnet"))})
#####################################
# EMPTY 
#####################################
#' @rdname empty
#' @aliases empty,mgnet
setMethod("empty", c("mgnet"),function(object){
  length(object@data)==0 & length(object@meta)==0 & length(object@taxa)==0 & 
    length(object@netw)==0 & length(object@comm)==0
})
################################################################################
################################################################################
# END BASE METHODS
################################################################################




################################################################################
################################################################################
# VERTICES SELECTION
################################################################################
################################################################################
#' Selection a subset of taxa from mgnet object.
#' 
#' @description The function takes as input a vector of logical or position 
#' indices to evaluate which taxa to keep. The function make a subset of the mgnet
#' object cutting the columns of data and the rows of taxa, making a subset of 
#' the network and creating a new communities object with only the conserved taxa.
#' 
#' @usage selection_vertices(object, idx)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param idx (Required) Vector of integer position indices or logical or names 
#' of taxa (like to \code{[} extractor function).
#'
#' @importFrom igraph subgraph
#' 
#' @rdname selection_vertices-methods
#' @docType methods
#' @export
setGeneric("selection_vertices", function(object,idx) standardGeneric("selection_vertices"))
#' @rdname selection_vertices-methods
#' @aliases selection_vertices,mgnet,vector
setMethod("selection_vertices", c("mgnet","vector"),
          function(object,idx){
            
            if(length(idx)>ntaxa(object)) stop("there are more indices than taxa")
            
            ifelse(length(object@data)!=0, data.new<-object@data[,idx,drop=F],  data.new<-object@data)
            ifelse(length(object@taxa)!=0, taxa.new<-object@taxa[idx,,drop=F],  taxa.new<-object@taxa)
            
            if(length(object@netw)!=0){
              netw.new<-igraph::subgraph(object@netw,idx)
            }else{
                netw.new<-object@netw
            }
              
            if(length(object@comm)!=0){
              comm.new <- object@comm
              if(is.character(idx)) idx <- which(taxaID(object)%in%idx)
              comm.new$membership <- object@comm$membership[idx]
              comm.new$vcount <- length(comm.new$membership)
              comm.new$modularity <- NA
            } else {
              comm.new <- object@comm
            }
            
            return(mgnet(data=data.new,
                         meta=object@meta,
                         taxa=taxa.new,
                         netw=netw.new,
                         comm=comm.new))
          })
#' @rdname selection_vertices-methods
#' @aliases selection_vertices,list,vector
setMethod("selection_vertices",c("list","vector"),
          function(object,idx){
            lapply(object, selectMethod(f="selection_vertices",signature=c("mgnet","vector")),
                   idx=idx)}
)
################################################################################
################################################################################
# END VERTICES SELECTION
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
#'  
#'@rdname arrange_vertices-methods
#'@docType methods
#'@export
setGeneric(name="arrange_vertices",
           def=function(obj, new_taxa) standardGeneric("arrange_vertices"))
#' @rdname arrange_vertices-methods
#' @aliases arrange_vertices,mgnet,matrix
setMethod("arrange_vertices",c("mgnet","matrix"),
          function(obj, new_taxa){
            
            if(length(obj@data)==0 | length(obj@taxa)==0 | length(obj@netw)==0) stop("data, taxa and netw cannot be empty")
            if(any(!(taxaID(obj)%in%rownames(new_taxa)))) stop("find at least a taxa not present in new_taxa")
            if(any(ranks(obj)!=colnames(new_taxa))) stop('new_taxa must have the same ranks as obj')
            
            sample_name <- sample_name(obj)
            ntaxa <- nrow(new_taxa)
            taxa_name <- rownames(new_taxa)
            nsample <- nsample(obj)
            
            data <- matrix(0,nrow=nsample,ncol=ntaxa,
                           dimnames=list(sample_name,taxa_name))
            data[,taxaID(obj)] <- obj@data
            
            adj <- as_adjacency_matrix(obj@netw,sparse=F,attr="weight")
            adj.new <- matrix(0,nrow=ntaxa,ncol=ntaxa,
                              dimnames=list(taxa_name,taxa_name))
            adj.new[taxaID(obj),taxaID(obj)] <- adj
            
            return(mgnet(data=data,taxa=new_taxa,adj=adj.new))
          })
#' @rdname arrange_vertices-methods
#' @aliases arrange_vertices,list,matrix
setMethod("arrange_vertices",c("list","matrix"),
          function(obj,new_taxa){
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
#' @description allows you to remove communities based on the number of 
#' vertices.
#' 
#' @param obj \code{\link{mgnet-class}}
#' @param size integer indicates the vertex number threshold
#' @param trim logical. If true, the function removes all nodes not belonging 
#' to a community with size equal to size. If false the filtered vertices are 
#' set as isolated.
#'  
#' @importFrom igraph vcount induced.subgraph V
#' @rdname remove_smaller_comm-methods
#' @docType methods
#' @export
setGeneric("remove_smaller_comm", function(obj,size,trim) standardGeneric("remove_smaller_comm"))
#' @rdname remove_smaller_comm-methods
#' @aliases remove_smaller_comm,mgnet,logical
setMethod("remove_smaller_comm", c("mgnet","numeric","logical"), function(obj, size, trim){
  
  #Checks Arguments
  if(!is.numeric(size) | round(size)!=size | size<0) stop("size must be integer greater than 0")
  if(!is.logical(trim)) stop("keep must be logical")
  #End Checks
  
  graph <- obj@netw
  comm <- obj@comm
  
  keep.comm.names <- as.numeric(names(sizes(comm)[sizes(comm)>=size]))
  keep.comm.vids <-  which(comm$membership %in% keep.comm.names)
  
  if(!trim){
    data <- obj@data
    taxa <- obj@taxa
    graph.sub <- graph
    comm.sub <- comm
    comm.sub$membership[setdiff(1:comm$vcount, keep.comm.vids)] <- 0
    comm.sub$modularity <- NA
  } else {
    data <- obj@data[,keep.comm.vids]
    taxa <- obj@taxa[keep.comm.vids,]
    graph.sub <- induced.subgraph(graph, V(graph)[keep.comm.vids])
    comm.sub <- comm
    comm.sub$membership <- comm$membership[keep.comm.vids]
    comm.sub$vcount <- length(comm.sub$membership)
    comm.sub$modularity <- NA
  }
  
  
  return(new("mgnet",
             data=data, meta=obj@meta, taxa=taxa,
             netw=graph.sub,comm=comm.sub))
})
#' @rdname remove_smaller_comm-methods
#' @aliases remove_smaller_comm-methods,list,numeric,logical
setMethod("remove_smaller_comm",c("list","numeric","logical"),
          function(obj,size,trim){
            lapply(obj, selectMethod(f="remove_smaller_comm",
                                     signature=c("mgnet","numeric","logical")),
                   size=size,trim=trim)})
################################################################################
################################################################################
# END REMOVE SMALLER COMMUNITIES
################################################################################
################################################################################




################################################################################
################################################################################
# GRAPHICAL DECORATIONS
################################################################################
################################################################################
#' Add graphical decorations
#' 
#'@description Add graphical descriptions to network:
#'\itemize{
#'  \item Vertex size proportional to mean clr abundances over samples.
#'  \item Link color could be red (+) or blue (-) respect the weights sign.
#'}
#' 
#' @param obj \code{\link{mgnet-class}}
#'  
#' @rdname default_decoration-methods
#' @docType methods
#' @export
setGeneric("default_decoration", function(obj) standardGeneric("default_decoration"))
#' @importFrom igraph V<- E<-
#' @importFrom grDevices rgb
#' @rdname default_decoration-methods
#' @aliases default_decoration,mgnet
setMethod("default_decoration", "mgnet", function(obj){
  
  if(length(obj@data)==0 | length(obj@netw)==0) stop("slot data and netw must be present")
  
  # Vertex size
  V(obj@netw)$size  <- 4 + colMeans(clr(obj@data+1)) + abs(min(colMeans(clr(obj@data+1))))
  
  # Edges color and width
  w <- E(obj@netw)$weight
  E(obj@netw)$color <- ifelse(w>0, grDevices::rgb(0,0,1,.5), grDevices::rgb(1,0,0,.5))
  E(obj@netw)$width <- abs(w) / max(abs(w))
  
  return(obj)
})
#' @rdname default_decoration-methods
#' @aliases default_decoration-methods,list,numeric,logical
setMethod("default_decoration","list",
          function(obj){
            lapply(obj, selectMethod(f="default_decoration",
                                     signature="mgnet"))})
################################################################################
################################################################################
# GRAPHICAL DECORATIONS
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
#' @rdname default_decoration-methods
#' @aliases default_decoration,igraph
setMethod("layout_signed","igraph",function(obj){
  
  graph.sub <- subgraph.edges(graph=obj,
                              eids=which(E(obj)$weight>0),
                              delete.vertices=FALSE)
  
  layout <- layout.fruchterman.reingold(graph.sub)
  return(layout)
})
#' @rdname default_decoration-methods
#' @aliases default_decoration,mgnet
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
#'@export
plot.mgnet <- function(x,layout,...) {
  if(missing(layout)){
    layout <- layout_signed(netw(x))
  }
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
#' @param sign character indicates the sign of the edges. Possible values are
#' "positive","negative","all".
#' @param type character with possible values "intra","extra","all".
#'
#' @importFrom stats setNames
#' @importFrom igraph degree subgraph.edges intersection crossing
#' @rdname degree_mgnet-methods
#' @docType methods
#' @export
setGeneric("degree_mgnet", function(obj,sign,type) standardGeneric("degree_mgnet"))
#' @rdname degree_mgnet-methods
#' @aliases degree_mgnet,mgnet,character,character
setMethod("degree_mgnet",c("mgnet","character","character"),function(obj,sign,type){
 
  if(length(comm(obj))==0) stop("communities slots is empty")
  
  sign <- match.arg(sign,c("positive","negative","all"))
  type <- match.arg(type,c("intra","extra","all"))
  
  
  n <- netw(obj)
  c <- comm(obj)
  
  # Modify isolated nodes as igraph want
  #------------------------------------#
  c$membership[c$membership==0] <- (length(c)+1):((length(c)+1)+length(c$membership[c$membership==0])-1)
  
  if(sign=="positive"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight>0],delete.vertices=FALSE)
  } else if(type=="extra"){
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
#' @aliases degree_mgnet-methods,list,character,character
setMethod("degree_mgnet",c("list","character","character"),
          function(obj,sign,type){
            lapply(obj, selectMethod(f="degree_mgnet",
                                     signature=c("mgnet","character","character")),
                   sign=sign,type=type)})
#' @rdname degree_mgnet-methods
#' @aliases degree_mgnet,mgnet,character,missing
setMethod("degree_mgnet",c("mgnet","character","missing"),function(obj,sign,type){
  
  sign <- match.arg(sign,c("positive","negative","all"))

  n <- netw(obj)
  c <- comm(obj)

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
  
  return(degree(sub.sign))
})
#' @rdname degree_mgnet-methods
#' @aliases degree_mgnet-methods,list,character,missing
setMethod("degree_mgnet",c("list","character","missing"),
          function(obj,sign,type){
            lapply(obj, selectMethod(f="degree_mgnet",
                                     signature=c("mgnet","character","missing")),
                   sign=sign)})
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
#' @param sign character indicates the sign of the edges. Possible values are
#' "positive","negative","all".
#' @param type character with possible values "intra","extra","all".
#'
#' @importFrom igraph degree subgraph.edges intersection
#' @rdname strength_mgnet-methods
#' @docType methods
#' @export
setGeneric("strength_mgnet", function(obj,sign,type) standardGeneric("strength_mgnet"))
#' @rdname strength_mgnet-methods
#' @aliases strength_mgnet,mgnet,character,character
setMethod("strength_mgnet",c("mgnet","character","character"),function(obj,sign,type){
  
  if(length(comm(obj))==0) stop("communities slots is empty")
  
  sign <- match.arg(sign,c("positive","negative","all"))
  type <- match.arg(type,c("intra","extra","all"))
  
  n <- netw(obj)
  c <- comm(obj)
  
  # Modify isolated nodes as igraph want
  #------------------------------------#
  c$membership[c$membership==0] <- (length(c)+1):((length(c)+1)+length(c$membership[c$membership==0])-1)
  
  if(sign=="positive"){
    sub.sign <- subgraph.edges(graph=n,eids=E(n)[E(n)$weight>0],delete.vertices=FALSE)
  } else if(type=="extra"){
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
#' @aliases strength_mgnet-methods,list,character,character
setMethod("strength_mgnet",c("list","character","character"),
          function(obj,sign,type){
            lapply(obj, selectMethod(f="strength_mgnet",
                                     signature=c("mgnet","character","character")),
                   sign=sign,type=type)})
#' @rdname strength_mgnet-methods
#' @aliases strength_mgnet,mgnet,character,missing
setMethod("strength_mgnet",c("mgnet","character","missing"),function(obj,sign,type){
  
  sign <- match.arg(sign,c("positive","negative","all"))
  n <- netw(obj)
  c <- comm(obj)

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
  
  adj.sub <- as_adjacency_matrix(sub.sign,attr="weight",sparse=FALSE)
  strength <- setNames(rowSums(adj.sub),taxaID(obj))
  return(strength)
})
#' @rdname strength_mgnet-methods
#' @aliases strength_mgnet-methods,list,character,missing
setMethod("strength_mgnet",c("list","character","missing"),
          function(obj,sign,type){
            lapply(obj, selectMethod(f="strength_mgnet",
                                     signature=c("mgnet","character","missing")),
                   sign=sign)})
################################################################################
################################################################################
# END STRENGTH MGNET
################################################################################
################################################################################