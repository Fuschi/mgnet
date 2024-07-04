#' Calculate Signed and Community-Specific Vertex Degrees in mgnet Objects
#'
#' This function computes the degrees of vertices in an `mgnet` object based on the signs of edge weights
#' and the specified relationship to community structure. It allows for focusing on positive, negative,
#' or all edges, and can differentiate between edges within the same community (intra-community) and
#' edges between different communities (inter-community).
#'
#' @param obj An `mgnet` object containing a network structured as an `igraph` object.
#' @param sign A character string indicating the sign of the edges to consider.
#'             Accepted values are "positive", "negative", or "all".
#' @param type A character string indicating the type of community relationship to consider for the
#'             degree calculation. Accepted values are "intra" for edges within the same community,
#'             "extra" for edges between different communities, or "all" for no community-based filtering.
#'
#' @return An integer vector containing the degrees of the vertices. Each element in the vector
#'         corresponds to a vertex in the `mgnet` object, ordered according to the vertex order in the
#'         igraph object.
#'
#' @details
#' The degree of a vertex is the number of edges connected to it, with possible restrictions based on
#' edge sign and community membership:
#' - `sign = "positive"`: Only edges with positive weights are considered.
#' - `sign = "negative"`: Only edges with negative weights are considered.
#' - `type = "intra"`: Only edges within the same community as the vertex are considered.
#' - `type = "extra"`: Only edges between different communities are considered.
#'
#' This function is particularly useful for analyzing the connectivity patterns in ecological networks
#' or social networks where edges are labeled with positive or negative weights to represent different
#' types of interactions or relationships, and where community structure is an important aspect.
#'
#' @importFrom igraph degree delete_edges E crossing
#' @export
#' @seealso \code{\link[igraph]{degree}} for the base function used to compute vertex degrees.
setGeneric("degree_mgnet", function(obj, sign = "all", type = "all") standardGeneric("degree_mgnet"))

setMethod("degree_mgnet", "mgnet",
          function(obj, sign = "all", type = "all") {
            sign <- match.arg(sign, c("positive", "negative", "all"))
            type <- match.arg(type, c("intra", "extra", "all"))
            
            # Ensure the network is present
            if (is.null(obj@network)) {
              stop("The mgnet object does not contain a network.")
            }
            
            # Ensure community data is available if needed
            if (type != "all" && is.null(obj@community)) {
              stop("Community data is required for intra or extra community calculations.")
            }
            
            net <- obj@network
            com <- as_igraph_communities(obj@community)
            
            # Filter edges based on sign (links to remove)
            eids_sign <- E(net)$weight * switch(sign, positive = 1, negative = -1, all = 0) < 0

            # Filter edges based on community intra/extra types (links to remove)
            eids_type <- switch(type, 
                                intra = igraph::crossing(com, net), 
                                extra = !igraph::crossing(com, net), 
                                all = rep(FALSE, vcount(net))) 

            # Remove links
            sub_netw <- delete_edges(net, which(eids_sign | eids_type) )
            
            # Calculate degree
            return(igraph::degree(sub_netw))
          })

setMethod("degree_mgnet", "mgnetList",
          function(obj, sign = "all", type = "all") {
            lapply(obj, function(x) degree_mgnet(x, sign = sign, type = type))
          })






















#' Signed Vertices Degrees with Communities Info
#'
#' Calculates the vertices degree using the information on the sign of the edges and the communities.
#' 
#' @param obj An mgnet object.
#' @param sign Character indicating the sign of the edges to consider. 
#' Possible values: "positive" for positive weights, "negative" for negative weights, or "all" for all weights.
#' @param type Character indicating the type of edges to include. 
#' Possible values: "intra" for edges within communities, "extra" for edges between communities, or "all" for all edges.
#' @return A vector of vertex degrees.
#' @importFrom igraph crossing intersection degree
#' @export
#' @aliases degree_mgnet,mgnet-method degree_mgnet,mgnetList-method
setGeneric("degree_mgnet", function(obj, sign = "all", type = "all") standardGeneric("degree_mgnet"))

setMethod("degree_mgnet", "mgnet",
          function(obj, sign = "all", type = "all") {
            sign <- match.arg(sign, c("positive", "negative", "all"))
            type <- match.arg(type, c("intra", "extra", "all"))
            
            if (type != "all" && length(community(obj)) == 0) {
              stop("With 'type' set, the 'community(' slot must not be empty.")
            }
            
            network <- network(obj)
            communities <- community(obj)
            
            if (sign == "positive") {
              subgraph <- subgraph.edges(graph = network, eids = E(network)[E(network)$weight > 0], delete.vertices = FALSE)
            } else if (sign == "negative") {
              subgraph <- subgraph.edges(graph = network, eids = E(network)[E(network)$weight < 0], delete.vertices = FALSE)
            } else {
              subgraph <- network
            }
            
            if (type != "all") {
              communities$membership[communities$membership == 0] <- seq_along(communities$membership)[communities$membership == 0] + max(communities$membership)
            }
            
            if (type == "intra") {
              subgraph <- subgraph.edges(graph = network, eids = E(network)[!crossing(communities, network)], delete.vertices = FALSE)
            } else if (type == "extra") {
              subgraph <- subgraph.edges(graph = network, eids = E(network)[crossing(communities, network)], delete.vertices = FALSE)
            }
            
            intersection_graph <- intersection(subgraph, network)
            return(degree(intersection_graph))
          })

setMethod("degree_mgnet", "mgnetList",
          function(obj, sign = "all", type = "all") {
            lapply(obj, function(x) degree_mgnet(x, sign = sign, type = type))
          })