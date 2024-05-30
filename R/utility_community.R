#' Convert igraph Communities to mgnet Format
#'
#' Converts an `igraph` communities object to the `mgnet` format by reassigning isolated nodes
#' to a community labeled as '0'. Unlike `igraph`, which treats isolated nodes as individual communities,
#' `mgnet` groups all isolated nodes into a single community for analysis purposes.
#'
#' @param comm An object of class `communities`, typically obtained from community detection algorithms in `igraph`.
#'
#' @return The modified `comm` object with isolated nodes reassigned to community '0'.
#'         The structure of `comm` is preserved, with changes applied only to the membership vector.
#'
#' @export
#' @importFrom igraph communities sizes
as_mgnet_communities <- function(comm) {
  if (!inherits(comm, "communities")) {
    stop("comm must belong to communities class", call. = FALSE)
  }
  
  # Identify isolated nodes based on community sizes
  isolated_memberships <- which(sizes(comm) == 1)
  # Reassign membership of isolated nodes to 0
  comm$membership[comm$membership %in% isolated_memberships] <- 0
  
  return(comm)
}


#' Convert mgnet Communities to igraph Format
#'
#' Transforms a communities object from the `mgnet` format back to `igraph` by assigning
#' unique community numbers to isolated nodes, which are initially grouped under '0' in `mgnet`.
#' This operation facilitates compatibility with `igraph`'s community analysis functions by ensuring
#' each isolated node is treated as its own community.
#'
#' @param comm An object of class `communities`, typically after manipulation or analysis within an `mgnet` context.
#'
#' @return The modified `comm` object with isolated nodes assigned unique community numbers, 
#'         starting from the highest existing community number + 1. The overall structure of `comm`
#'         remains unchanged, with modifications applied to the membership vector to reflect `igraph` conventions.
#'
#' @export
#' @importFrom igraph communities membership
as_igraph_communities <- function(comm) {
  if (!inherits(comm, "communities")) {
    stop("comm must belong to communities class", call. = FALSE)
  }
  
  isolated <- which(membership(comm) == 0)
  
  # If there are isolated nodes to reassign
  if (length(isolated) > 0) {
    new_membership <- comm$membership
    max_community <- max(new_membership)
    # Generate new community numbers for isolated nodes
    new_communities <- seq(from = max_community + 1, 
                           by = 1, length.out = length(isolated))
    new_membership[isolated] <- new_communities
    
    comm$membership <- new_membership
  }
  
  return(comm)
}