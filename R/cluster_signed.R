#' Get Communities of a Signed Weighted Graph
#'
#'@description Adaptation in R of the algorithm for detention of weighted 
#' communities with sign of a graph developed by Sergio Gomez.
#' https://deim.urv.cat/~sergio.gomez/radatools.php. 
#' The algorithm tries to find dense subgraph, also called communities via
#' optimization of a signed definition of modularity score.
#' 
#'@param obj weighted undirected network belong to \code{\link{igraph}} class 
#'or an mgnet.
#'@param OS (default Linux) string with the operating system running. Possible choices are 
#'"Linux","Windows","Mac".
#'@param Resistance resistance of nodes to join communities, as a common 
#'self-loop positive or negative real number, default no resistance with value set
#'to 0.
#'@param Penalty_Coefficien relative importance of null-case term non-negative 
#'real number default set to 1.
#'@param add.names logical with default value set to TRUE. It indicates whether 
#'you want to use the name of the vertices also in the resulting communities or 
#'simply numerical indexes.
#'
#'@return \code{\link{communities}} igraph object able to manage to communities
#'graph info. The unique difference from igraph routine is the treatment with 
#'the isolated nodes. In this case all isolated nodes are classified in the
#'community \code{'0'} and not as different communities of size one.
#'
#'@importFrom igraph V graph_from_adjacency_matrix write_graph make_clusters is.weighted is.directed as_adjacency_matrix make_clusters sizes is_named
#'@importFrom stringr str_split
#'
#' @rdname cluster_signed-methods
#' @docType methods
#' @export
setGeneric("cluster_signed", function(obj,OS="Linux", 
                                      Resistance=0, Penalty_Coefficien=1,
                                      add.names=TRUE) standardGeneric("cluster_signed"))
#' @rdname cluster_signed-methods
setMethod("cluster_signed", "igraph", function(obj,OS="Linux",
                                               Resistance=0, Penalty_Coefficien=1,
                                               add.names=TRUE){
  
  graph <- obj
  #Check Graph
  OS <- match.arg(OS,c("Linux","Windows","Mac"))
  if(!is.weighted(graph)) stop("graph must be weighted graph")
  if(is.directed(graph))  stop("graph must be undirected")
  if(file.exists("graph.net")) file.remove("graph.net")
  if(!is.numeric(Resistance)) stop("Resistance must be numeric")
  if(!is.numeric(Penalty_Coefficien)) stop("Penalty_Coefficien must be a number >= 0")
  if(Penalty_Coefficien<0) stop("Penalty_Coefficien must be a number >= 0")
  if(add.names & !igraph::is_named(obj)) stop("graph has not vertices names attribute")
  
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
  cmd <- paste(cmd, Resistance, Penalty_Coefficien, path.graph, path.result)
  
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
#' @rdname cluster_signed-methods
setMethod("cluster_signed","mgnet",
          function(obj,OS="Linux",add.names=TRUE){
            comm <- cluster_signed(netw(obj),OS=OS,add.names=add.names)
            return(mgnet(data=obj@data, meta_sample=obj@meta_sample,
                         taxa=obj@taxa, meta_taxa=obj@meta_taxa,
                         log_data=obj@log_data,
                         netw=obj@netw, comm=comm))
          })
#' @rdname cluster_signed-methods
setMethod("cluster_signed","list",
          function(obj,OS="Linux",add.names=TRUE){
            return(lapply(obj, selectMethod(f="cluster_signed",signature="mgnet"),
                   OS=OS,add.names=add.names))})
