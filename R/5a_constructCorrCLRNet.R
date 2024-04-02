# CONSTRUCT CORRELATION CLR NETWORK
#------------------------------------------------------------------------------#
#' Construct Correlation Network from CLR-Transformed Data
#'
#' Constructs a correlation network from compositional data using centered log-ratio (CLR) 
#' transformation followed by correlation analysis. This function is designed to handle 
#' compositional data's inherent constraints by applying CLR transformation, facilitating 
#' meaningful correlation analysis. Users can choose from various zero replacement strategies, 
#' CLR methods, and thresholding techniques to tailor the network construction process to their 
#' specific data and research questions.
#'
#' @param object An `mgnet` or `mgnetList` object containing abundance data.
#' @param zero_strategy Strategy for handling zeros in the abundance data before CLR transformation. 
#' Options are "unif" for replacing zeros with values from a uniform distribution and "const" 
#' for a constant value replacement. Default is "unif".
#' @param clr_method Specifies the CLR transformation method: "clr" for standard CLR, "iclr" for 
#' interquartile CLR, or "none" to bypass CLR transformation and use precomputed norm_abundance data. Default is "clr".
#' @param cor_method Correlation method ("pearson", "spearman", "kendall") used for computing 
#' pairwise correlations between taxa after CLR transformation. Default is "pearson".
#' @param thresh_method Method for thresholding the correlation matrix to construct the adjacency 
#' matrix: "absolute" applies an absolute value threshold, "density" targets a specific edge 
#' density, and "p-value" uses p-value adjustment to threshold correlations. 
#' @param thresh_value Numeric threshold applied according to `thresh_method`. For "absolute" and 
#' "p-value", this represents the minimum correlation or adjusted p-value for inclusion in the 
#' network. For "density", this represents the desired edge density. 
#' @param padj_method Method for adjusting p-values when `thresh_method` is "p-value". Options 
#' include "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), "BY" 
#' (Benjamini-Yekutieli), "fdr" (false discovery rate), and "none" (no adjustment). Default is "none".
#'
#' @return The input `mgnet` or `mgnetList` object with an updated `network` slot containing the 
#' constructed correlation network. The network is represented as an `igraph` object for `mgnet` 
#' and a list of `igraph` objects for `mgnetList`.
#'
#' @details
#' The function operates in several steps to transform the abundance data into a correlation 
#' network:
#' 1. **Zero Replacement:** Addressing the issue of zeros in compositional data, which can skew 
#' CLR transformation and correlation analysis. The chosen `zero_strategy` is applied across 
#' all samples.
#' 2. **CLR Transformation:** Normalizing the data using CLR transformation (`clr`) or 
#' interquartile CLR (`iclr`), making it suitable for correlation analysis by removing 
#' compositional constraints.
#' 3. **Correlation Analysis:** Computing pairwise correlations between taxa using the specified 
#' `cor_method`, resulting in a correlation matrix.
#' 4. **Threshold Application:** Constructing the adjacency matrix by applying the selected 
#' threshold method (`thresh_method`) and value (`thresh_value`). This step determines which 
#' correlations are considered significant enough to constitute edges in the network.
#' 5. **Network Construction:** The adjacency matrix is used to construct a weighted, undirected 
#' correlation network, represented as an `igraph` object.
#'
#' This approach enables the construction of meaningful correlation networks from compositional 
#' data, facilitating the exploration of complex microbial community interactions and structure.
#'
#'
#' @seealso
#' \code{\link[mgnet]{clr}}, \code{\link[mgnet]{iclr}} for centered-log ratio details.
#' \code{\link[mgnet]{zero_dealing}} for details on zero replacement strategies.
#'
#' @export
#' @importFrom igraph graph_from_adjacency_matrix make_empty_graph
#' @importFrom methods validObject
#' @name constructCorrCLRNet
#' 
#' @aliases constructCorrCLRNet,mgnet-method constructCorrCLRNet,mgnetList-method
#' @export
setGeneric("constructCorrCLRNet", function(object, zero_strategy = "unif", clr_method ="clr",
                                           cor_method = "pearson", thresh_method, thresh_value, padj_method="none") standardGeneric("constructCorrCLRNet"))

setMethod("constructCorrCLRNet", "mgnet",
          function(object, zero_strategy = "unif", clr_method ="clr",
                   cor_method = "pearson", thresh_method, thresh_value, padj_method="none"){
            
            
            # Checks
            zero_strategy <- match.arg(zero_strategy,c("unif","const"))
            clr_method <- match.arg(clr_method,c("clr","iclr","none"))
            cor_method <- match.arg(cor_method,c("pearson","spearman","kendall"))
            thresh_method <- match.arg(thresh_method,c("absolute","density","p-value"))
            padj_method <- match.arg(padj_method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
            if(thresh_method != "p-value" && padj_method != "none") stop("adjust different from 'none' can be set only for thresh_method equal to p-value")
            if(!is.numeric(thresh_value) | thresh_value<0 | thresh_value>1) stop("thresh_value must number in range [0,1]")
            if( clr_method != "none" && length(object@norm_abundance)!=0 ) stop("norm_abundance matrix is not empty. Please set clr_method to 'none' or remove from the object.")
            if( clr_method != "none" && length(object@abundance)==0 ) stop("abundance matrix missing and i cannot calculate the clr") 

            if( clr_method != "none" && length(object@norm_abundance)==0 ){
              
              norm_abundance <- if(clr_method == "clr") {
                  mgnet::clr(zero_dealing(object@abundance, method = zero_strategy))
                } else {
                  mgnet::iclr(zero_dealing(object@abundance, method = zero_strategy))
                }
              
            } else if ( clr_method == "none" && length(object@norm_abundance)!=0 ){
              
              norm_abundance <- object@norm_abundance
              
            } else {
              
              stop("why are you here?")
              
            }
            
            if(thresh_method=="absolute"){
              
              adj <- cor(norm_abundance, method=cor_method)
              adj <- adj * (abs(adj) >= thresh_value)
              diag(adj) <- 0
              
            } else if(thresh_method=="density"){
              
              adj <- mgnet::adjacency_edge_density(x = norm_abundance, method = cor_method, th = thresh_value)
              
            } else if(thresh_method=="p-value"){
              
              adj <- mgnet::adjacency_p_adjust(x = norm_abundance, method = cor_method, adjust = padj_method, alpha = thresh_value)
              
            }
            
            diag(adj) <- 0
            if(all(adj==0)){
              network <- make_empty_graph( n = ntaxa(object), directed = FALSE)
              V(network)$name <- taxa_id(object)
              E(network)$weight <- numeric(0)
            } else {
              network <- igraph::graph_from_adjacency_matrix(adjmatrix = adj, 
                                                             mode = "undirected", 
                                                             weighted = TRUE, 
                                                             diag = FALSE)
            }
            
            object@norm_abundance <- norm_abundance
            object@network <- network
            validObject(object)
            return(object)
          })


setMethod("constructCorrCLRNet", "mgnetList",
          function(object, zero_strategy = "unif", clr_method ="clr",
                   cor_method = "pearson", thresh_method, thresh_value, padj_method="none"){
            
            object@mgnets <- sapply(object@mgnets, function(x){
              constructCorrCLRNet(x,
                                  zero_strategy = "unif", clr_method ="clr",
                                  cor_method = "pearson", thresh_method, thresh_value, padj_method="none")},
              simplify = FALSE, USE.NAMES = TRUE)
            
            validObject(object)
            return(object)
          })