#' Community Detection in Signed Weighted Graphs
#'
#' @description
#' Implements a wrapper for community detection in signed weighted graphs,
#' adapting the algorithm developed by Sergio Gomez. This method is designed
#' to identify dense subgraphs, known as communities, by optimizing a signed
#' modularity score. The algorithm accommodates weighted undirected graphs,
#' applying a threshold to identify the most significant correlations as edges.
#'
#' @param object An `igraph` object representing a weighted undirected network or
#' an `mgnet` object containing such a network or an `mgnetList`.
#' @param resistance (Default: 0) Node resistance to joining communities,
#' specified as a positive or negative real number. This parameter simulates
#' a common self-loop resistance across nodes, with the default indicating no
#' resistance.
#' @param penalty (Default: 1) Specifies the relative importance of
#' the null-case term as a non-negative real number. This coefficient affects
#' how strongly the algorithm penalizes or favors the formation of communities.
#' @param add.names (Default: FALSE) Determines whether vertex names (if available)
#' are used in the output, associating each node with its community assignment.
#' When FALSE, community assignments are returned as numerical indices only.
#'
#' @return An object of class `communities` as defined in the `igraph` package.
#' This object contains community assignments for each node in the input graph.
#' Isolated nodes are grouped into a special community labeled '0', rather than
#' being considered as separate communities of size one.
#'
#' @details
#' The function is an adaptation of Sergio Gomez's Radatools community detection
#' algorithm for signed networks. It identifies community structures by optimizing
#' a signed version of the modularity score, accounting for the positive and
#' negative weights of edges. The algorithm is particularly suited for analyzing
#' networks where relationships can have both positive (affiliative) and negative
#' (antagonistic) connotations.
#'
#' @references
#' Gomez, S. et al. (2009). Analysis of community structure in networks of correlated data. Phys. Rev. E. 
#'
#' @seealso
#' \code{\link[igraph]{communities}}, for information on handling community objects in `igraph`.
#'
#' @export
#' @importFrom igraph graph_from_adjacency_matrix write_graph make_clusters is.weighted is.directed as_adjacency_matrix make_clusters sizes is_named
#' @importFrom stringr str_split
#' @name cluster_signed
#' @aliases cluster_signed,igraph-method cluster_signed,mgnet-method cluster_signed,mgnetList-method
setGeneric("cluster_signed", 
           function(object, resistance=0, penalty=1, add.names=FALSE) standardGeneric("cluster_signed"))


setMethod("cluster_signed", "igraph", 
          function(object,resistance=0, penalty=1, add.names=FALSE){
  
  # Automatically detect the operating system
  detectedOS <- Sys.info()["sysname"]
  if (detectedOS == "Windows") {
    OS <- "Windows"
  } else if (detectedOS == "Darwin") {  # Sys.info() returns 'Darwin' for macOS
    OS <- "Mac"
  } else if (detectedOS == "Windows") {
    # Assuming Linux for all non-Windows and non-Mac systems
    OS <- "Linux"
  } else {
    stop("What operating system are you using??")
  }
  
  graph <- object
  #Check Graph
  if(!is.weighted(graph) |! is.directed(graph)) stop("graph must be weighted unidrected graph")
  if(file.exists("graph.net")) file.remove("graph.net")
  if(!is.numeric(resistance)) stop("resistance must be numeric")
  if(!is.numeric(penalty)) stop("penalty must be a number >= 0")
  if(penalty<0) stop("penalty must be a number >= 0")
  if(add.names & !igraph::is_named(object)) stop("graph has not vertices names attribute")
  
  adj <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)
  #End Checks
  
  #get path of executable Communities_Detection.exe
  path <- system.file("exec", package="mgnet", mustWork=TRUE)
  path <- paste(path,"/",sep="")
  
  #write graph in pajek format as request from executable
  write_graph(graph  = graph_from_adjacency_matrix(adjmatrix=adj, mode="undirected", weighted=TRUE),
              file   = paste(path,"graph.net",sep=""),
              format ="pajek")
  
  #path for graph and results files
  path.graph <- paste(path,"graph.net",sep="")
  path.result <- paste(path, "res.txt", sep="")
  
  #save command for the execution
  if(OS=="Linux"){
    cmd <- paste(path,"Communities_Detection_Linux.exe none WS l 1",sep="")
  } else if(OS=="Windows"){
    cmd <- paste(path,"Communities_Detection_Windows.exe none WS l 1",sep="")
  } else if(OS=="Mac"){
    cmd <- paste(path,"Communities_Detection_Mac.exe none WS l 1",sep="")
  }
  cmd <- paste(cmd, resistance, penalty, path.graph, path.result)
  
  #make communities detection
  if(OS=="Windows"){
    system(cmd,ignore.stdout=TRUE)
  } else {
    system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
  }
  #open file of results
  res.file = file(path.result, "r")
  
  #read and store results
  file.lines <- readLines(res.file)
  info <- file.lines[2]
  modularity <- as.numeric(str_split(file.lines[3]," ")[[1]][3])
  
  vertex.num <- as.numeric(str_split(file.lines[5]," ")[[1]][4])
  comm.num <- as.numeric(str_split(file.lines[6]," ")[[1]][4])
  comm.vert <- vector(mode="list",length=comm.num)
  
  for(line in 8:(7+comm.num)){
    
    i <- line-7
    #comm[[i]]
    tmp <- unlist(str_split(file.lines[line]," "))
    tmp <- as.numeric(tmp[2:length(tmp)])
    
    comm.vert[[i]] <- tmp
  }
  
  comm <- vector(mode="integer", length=vertex.num)
  for(c in 1:comm.num){comm[comm.vert[[c]]] <- c}
  
  #close connection
  close(res.file)
  
  #remove files
  file.remove(path.graph)
  file.remove(path.result)
  
  res <- make_clusters(graph=graph,
                       membership=comm,
                       algorithm="signed weights louvain",
                       modularity=modularity)
  
  isolated.membership <- which(sizes(res)==1)
  res$membership[res$membership %in% isolated.membership] <- 0
  
  if(add.names) names(res$membership) <- V(graph)$name
  
  #return the results as communities structure of igraph package
  return(res)
})

setMethod("cluster_signed", "mgnet",
          function(object, resistance=0, penalty=1, add.names=TRUE){
            
            if(length(object@network)==0) stop("network missing")
            community(object) <- cluster_signed(network(object), 
                                               resistance=resistance, penalty=penalty, 
                                               add.names=add.names)
            validObject(object)
            return(object)
            
          })

setMethod("cluster_signed", "mgnetList",
          function(object, resistance=0, penalty=1, add.names=FALSE){
            
           object <- sapply(object, function(x) cluster_signed(x, 
                                                               resistance=resistance, penalty=penalty, 
                                                               add.names=add.names), 
                            simplify = FALSE, USE.NAMES = TRUE)
           validObject(object)
           return(object)
           })


