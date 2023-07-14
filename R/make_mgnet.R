################################################################################
################################################################################
# MAKE MGNET
################################################################################
################################################################################
#' Build a meta-genomic network.
#' 
#' @description The function has the purpose of building a meta-genomic network 
#' starting from a log_data slot in object following the following steps:
#' \itemize{
#' \item Calculates the pairwise taxa correlation matrix from the log_data slot.
#' \item Retrieves the adjacency matrix by applying a threshold to the lowest 
#' correlations with three possible criteria 'abosulte','density','p-value'.
#' \item Constructs the undirected, weighted, and signed network of class igraph.
#' }
#' 
#' @param mgnet mgnet object
#' @param cor.method correlation type, possible choices are "pearson","spearman","kendall".
#' @param thresh.method threshold method to retrieve the adjacency matrix. Possible
#' choices are "absolute","density","p-value".
#' @param adjust multiple test criteria, possible choices are "holm",
#' "hochberg","hommel","bonferroni","BH","BY","fdr","none".
#' @param thresh numeric indicates the threshold value.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' 
#' @rdname make_mgnet-methods
#' @docType methods
#' @export
setGeneric("make_mgnet", function(mgnet, cor.method, thresh.method, adjust, thresh) standardGeneric("make_mgnet"))
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("mgnet","character","character","character","numeric"),
          function(mgnet, cor.method, thresh.method, adjust, thresh){
            
            if(length(mgnet@log_data)==0) stop("data cannot be empty")
            
            # Check cor.method
            cor.method <- match.arg(cor.method,c("pearson","spearman","kendall"))
            # Check thresh.method
            thresh.method <- match.arg(thresh.method,c("absolute","density","p-value"))
            # Check adjust
            adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
            if(thresh.method!="p-value" & adjust!="none") stop("adjust different from none can be set only for thresh.method equal to p-value")
            # Check thresh
            if(!is.numeric(thresh) | thresh<0 | thresh>1) stop("thresh must number in range [0,1]")
            
            log_data <- mgnet@log_data
            
            if(thresh.method=="absolute"){
              adj <- cor(log_data,method=cor.method)
              adj <- adj * (abs(adj)>=thresh)
              diag(adj) <- 0
            } else if(thresh.method=="density"){
              adj <- mgnet::adjacency_edge_density(x=log_data,method=cor.method,th=thresh)
            } else if(thresh.method=="p-value"){
              adj <- mgnet::adjacency_p_adjust(x=log_data,method=cor.method,adjust=adjust,alpha=thresh)
            }
            
            return(mgnet(data=mgnet@data, meta_sample=mgnet@meta_sample,
                       taxa=mgnet@taxa, meta_taxa=mgnet@meta_taxa,
                       log_data=log_data, netw=adj))
          })
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("list","character","character","character","numeric"),
          function(mgnet, cor.method, thresh.method, adjust, thresh){
            lapply(mgnet, 
                   selectMethod(f="make_mgnet",signature=c("mgnet","character","character","character","numeric")),
                   cor.method=cor.method, thresh.method=thresh.method, adjust=adjust, thresh=thresh)
          })