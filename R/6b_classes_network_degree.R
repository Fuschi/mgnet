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