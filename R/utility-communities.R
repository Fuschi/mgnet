################################################################################
################################################################################
# AS MGNET COMMUNITIES
################################################################################
################################################################################
#' Convert igraph into mgnet communities.
#' 
#' @description The mgnet communities treat the isolated nodes differently respect 
#' igraph putting all of them in the community labeled as '0'. Contrary igraph
#' describes the isolated nodes as community of size 1. In conclusion the function
#' labels the isolate nodes as '0'.
#'
#' @param comm object belong to communities class.
#' 
#' @export
as_mgnet_communities <- function(comm){
  if(!is(comm,"communities")) stop("comm must belong to communities class")
  
  isolated.membership <- which(sizes(comm)==1)
  comm$membership[comm$membership %in% isolated.membership] <- 0
  return(comm)
}
################################################################################
################################################################################
# AS IGRAPH COMMUNITIES
################################################################################
################################################################################
#' Convert mgnet into igraph communities.
#' 
#' @description The mgnet communities treat the isolated nodes differently respect 
#' igraph putting all of them in the community labeled as '0'. Contrary igraph
#' describes the isolated nodes as community of size 1. In conclusion the function
#' modifies the the labels of the isolated nodes setting increasing number.
#'
#' @param comm object belong to communities class.
#' 
#' @export
as_igraph_communities <- function(comm){
  if(!is(comm,"communities")) stop("comm must belong to communities class")

  isolated <- which(membership(comm)==0)
  new.memb <- comm$membership
  ncomm <- names(sizes(comm))
  class(ncomm) <- "numeric"
  ncomm <- max(ncomm)
  new.memb[isolated] <- (ncomm+1):( length(new.memb[isolated])+ncomm )

  comm$membership <- new.memb
  return(comm)
}



