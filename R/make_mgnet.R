################################################################################
################################################################################
# MAKE MGNET
################################################################################
################################################################################
#' Build a meta-genomic network.
#' 
#' @description The function has the purpose of building a meta-genomic network 
#' starting from a mg object following the following steps:
#' \itemize{
#' \item calculates the clr abundances from data and if there are zeros, 
#' add the value 1 to all elements.
#' \item Evaluate the correlations between all possible pairs of taxa being able
#'  to choose between Pearson Kendall and Spearman.
#' \item derives the adjacency matrix by applying a threshold to the lowest 
#' correlations with three possible criteria. The simplest by giving it the 
#' absolute value, based on the density of the links or considering the p-value 
#' corrected with multiple tests.
#' \item constructs the undirected, weighted, and signed network of class igraph.
#' }
#' 
#' @param mgnet mgnet object
#' @param log.method log-ratio transformation choosen, possible choices are "base"
#' for the classical clr-transformation, "inter-quantile" for its interquantile version
#' and "stored" to use the log_data slot previous saved.
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
setGeneric("make_mgnet", function(mgnet, log.method, cor.method, 
                                  thresh.method, adjust, thresh) standardGeneric("make_mgnet"))
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("mgnet","character","character","character","character","numeric"),
          function(mgnet, log.method, cor.method, thresh.method, adjust, thresh){
            
            if(length(mgnet@data)==0) stop("data cannot be empty")

            # Check mgnet
            if(!isa(mgnet,"mgnet")) stop("mgnet must belong to mgnet class")
            # Check log.method
            log.method <- match.arg(log.method,c("base","inter-quantile","stored"))
            if(log.method=="stored" & length(mgnet@log_data)==0) stop("if log.method is stored the log_data slot cannot be empty")
            # Check cor.method
            cor.method <- match.arg(cor.method,c("pearson","spearman","kendall"))
            # Check thresh.method
            thresh.method <- match.arg(thresh.method,c("absolute","density","p-value"))
            # Check adjust
            adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
            if(thresh.method!="p-value" & adjust!="none") stop("adjust different from none can be set only for thresh.method equal to p-value")
            # Check thresh
            if(!is.numeric(thresh) | thresh<0 | thresh>1) stop("thresh must number in range [0,1]")
            
            if(log.method=="base"){
              mgnet <- save_log_data(mgnet,"base")
            } else if(log.method=="inter-quantile"){
              mgnet <- save_log_data(mgnet,"inter-quantile")
            } 
            
            log_data <- mgnet@log_data
            
            if(thresh.method=="absolute"){
              adj <- cor(log_data,method=cor.method)
              adj <- adj * (adj>=thresh)
            } else if(thresh.method=="density"){
              adj <- mgnet::adjacency_edge_density(x=log_data,method=cor.method,th=thresh)
            } else if(thresh.method=="p-value"){
              adj <- mgnet::adjacency_p_adjust(x=log_data,method=cor.method,adjust=adjust,alpha=thresh)
            }
            
            return(mgnet(data=mgnet@data, meta_sample=mgnet@meta_sample,
                         taxa=mgnet@taxa, meta_taxa=mgnet@meta_taxa,
                         log_data=log_data, adj=adj))
          })
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("list","character","character","character","character","numeric"),
          function(mgnet, log.method, cor.method, thresh.method, adjust, thresh){
            lapply(mgnet, 
                   selectMethod(f="make_mgnet",signature=c("mgnet","character","character","character","character","numeric")),
                   log.method=log.method, cor.method=cor.method, thresh.method=thresh.method,
                   adjust=adjust, thresh=thresh)
          })
