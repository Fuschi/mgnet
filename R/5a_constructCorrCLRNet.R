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
#' for a constant value replacement. Default is "const".
#' @param clr_method Specifies the CLR transformation method to be used. The default behavior
#' is to check for pre-computed normalized data stored in the `norm` slot of the `mgnet` object.
#' If normalized data is available, it uses this stored data unless another method is explicitly specified.
#' Available methods are:
#'   - \code{"stored"}: Use pre-computed normalized data from the `norm` slot. If no normalized data
#'     is present, an error is raised unless another method is specified.
#'   - \code{"clr"}: Perform standard CLR (Centered Log Ratio) transformation, suitable for most
#'     compositional data analysis scenarios. This will be used to calculate CLR from the raw abundance
#'     data if the `norm` slot is empty and no method is specified.
#'   - \code{"iclr"}: Perform an interquartile CLR transformation, which can be more robust to outliers.
#'     This method is selected only if explicitly specified.
#' The default setting automatically uses the `norm` slot if available; otherwise, it defaults to
#' calculating standard CLR data.
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
#' @param cores Number of cores to use for parallel processing; default is 1 (no parallel execution).
#'              This parameter is only applicable for `mgnetList` objects and is ignored for `mgnet` objects.
#'
#' @return The input `mgnet` or `mgnetList` object with an updated `netw` slot containing the
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
setGeneric("constructCorrCLRNet", function(object, zero_strategy = "const", clr_method = NULL,
                                           cor_method = "pearson", thresh_method = NULL,
                                           thresh_value = NULL, padj_method = NULL,
                                           cores = 1) standardGeneric("constructCorrCLRNet"))

setMethod("constructCorrCLRNet", "mgnet",
          function(object, zero_strategy = "const", clr_method = NULL,
                   cor_method = "pearson", thresh_method = NULL,
                   thresh_value = NULL, padj_method = NULL,
                   cores = 1){

            # Checks
            if(is.null(padj_method)) padj_method <- "none"
            if(is.null(clr_method) & length(norm(object))!=0) clr_method <- "stored"
            if(is.null(clr_method) & length(norm(object))==0) clr_method <- "clr"

            zero_strategy <- match.arg(zero_strategy,c("const","unif"))
            clr_method <- match.arg(clr_method,c("clr","iclr","stored"))
            cor_method <- match.arg(cor_method,c("pearson","spearman","kendall"))
            thresh_method <- match.arg(thresh_method,c("absolute","density","p-value"))
            padj_method <- match.arg(padj_method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
            if(thresh_method != "p-value" && padj_method != "none") stop("adjust different from 'none' can be set only for thresh_method equal to p-value")
            if(!is.numeric(thresh_value) | thresh_value<0 | thresh_value>1) stop("thresh_value must number in range [0,1]")
            if( clr_method == "stored" && length(object@norm)==0 ) stop("norm matrix missing")

            if( clr_method != "stored"){

              if(length(norm(object)) != 0) message("norm slot rewritten")

              norm <- if(clr_method == "clr") {
                  mgnet::clr(zero_dealing(object@abun, method = zero_strategy))
                } else {
                  mgnet::iclr(zero_dealing(object@abun, method = zero_strategy))
                }

            } else if ( clr_method == "stored" && length(object@norm)!=0 ){

              norm <- object@norm

            } else {

              stop("why are you here?")

            }

            if(thresh_method=="absolute"){

              adj <- cor(norm, method=cor_method)
              adj <- adj * (abs(adj) >= thresh_value)
              diag(adj) <- 0

            } else if(thresh_method=="density"){

              adj <- mgnet::adjacency_edge_density(x = norm, method = cor_method, th = thresh_value)

            } else if(thresh_method=="p-value"){

              adj <- mgnet::adjacency_p_adjust(x = norm, method = cor_method, adjust = padj_method, alpha = thresh_value)

            }

            diag(adj) <- 0
            if(all(adj==0)){
              netw <- make_empty_graph( n = ntaxa(object), directed = FALSE)
              V(netw)$name <- taxa_id(object)
              E(netw)$weight <- numeric(0)
            } else {
              netw <- igraph::graph_from_adjacency_matrix(adjmatrix = adj,
                                                             mode = "undirected",
                                                             weighted = TRUE,
                                                             diag = FALSE)
            }

            object@norm <- norm
            object@netw <- netw
            validObject(object)
            return(object)
          })


setMethod("constructCorrCLRNet", "mgnetList",
          function(object, zero_strategy = "const", clr_method = NULL,
                   cor_method = "pearson", thresh_method = NULL, thresh_value = NULL, padj_method = NULL,
                   cores = 1){
            
            if (cores > 1) {
              
              cl <- makeCluster(min(cores, detectCores()))
              on.exit(stopCluster(cl))

              object@mgnets <- parLapply(cl, object@mgnets, constructCorrCLRNet,
                                                   zero_strategy = zero_strategy, clr_method = clr_method,
                                                   cor_method = cor_method, thresh_method = thresh_method,
                                                   thresh_value = thresh_value, padj_method = padj_method)
            } else {
              object@mgnets <- lapply(object@mgnets, constructCorrCLRNet,
                                      zero_strategy = zero_strategy, clr_method = clr_method,
                                      cor_method = cor_method, thresh_method = thresh_method,
                                      thresh_value = thresh_value, padj_method = padj_method)
            }

            validObject(object)
            return(object)
          })