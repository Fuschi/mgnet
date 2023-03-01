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
<<<<<<< HEAD
#' @param mgnet mgnet object
=======
#' @param mg mg object
>>>>>>> main
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
<<<<<<< HEAD
setGeneric("make_mgnet", function(mgnet, cor.method, thresh.method, thresh, adjust=NULL) standardGeneric("make_mgnet"))
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("mgnet","character","character","numeric"),
          function(mgnet, cor.method, thresh.method, thresh, adjust=NULL){
            
            if(length(mgnet@data)==0) stop("data cannot be empty")
            if(!("geometric_mean"%in%sample_info(mgnet))) stop("sample geometric mean must be set before. See save_geometric_mean")
            
            # Check mgnet
            if(!isa(mgnet,"mgnet")) stop("mgnet must belong to mgnet class")
            # Check cor.method
            cor.method <- match.arg(cor.method,c("pearson","spearman","kendall"))
            # Check thresh.method
            thresh.method <- match.arg(thresh.method,c("absolute","density","p-value"))
            # Check adjust
            if(thresh.method!="p-value" & !is.null(adjust)) stop("adjust can be set only for thresh.method equal to p-value")
            if(thresh.method=="p-value"){
              adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
            }
            # Check thresh
            if(!is.numeric(thresh) | thresh<0 | thresh>1) stop("thresh must number in range [0,1]")
            
            
            if(thresh.method=="absolute"){
              adj <- cor(mgnet::CLR(mgnet),method=cor.method)
              adj <- adj * (adj>=thresh)
            } else if(thresh.method=="density"){
              adj <- mgnet::adjacency_edge_density(x=mgnet::CLR(mgnet),method=cor.method,th=thresh)
            } else if(thresh.method=="p-value"){
              adj <- mgnet::adjacency_p_adjust(x=mgnet::CLR(mgnet),method=cor.method,adjust=adjust,alpha=thresh)
            }
            
            return(mgnet(data=mgnet@data, meta_sample=mgnet@meta_sample,
                         taxa=mgnet@taxa, meta_taxa=mgnet@meta_taxa,
                         adj=adj))
          })
#' @rdname make_mgnet-methods
setMethod("make_mgnet", c("list","character","character","numeric"),
          function(mgnet, cor.method, thresh.method, thresh, adjust=NULL){
            lapply(mgnet, 
                   selectMethod(f="make_mgnet",signature=c("mgnet","character","character","numeric")),
                   cor.method=cor.method, thresh.method=thresh.method,
                   thresh=thresh, adjust=adjust)
          })
=======
setGeneric("make_mgnet", function(mg, cor.method, thresh.method, adjust, thresh) standardGeneric("make_mgnet"))
#' @rdname make_mgnet-methods
#' @aliases make_mgnet,mg,character,character,character,numeric
setMethod("make_mgnet", c("mg","character","character","character","numeric"),
          function(mg, cor.method, thresh.method, adjust, thresh){
            
  # Check mg
  if(!isa(mg,"mg")) stop("mg must belong to mg class")
  # Check cor.method
  cor.method <- match.arg(cor.method,c("pearson","spearman","kendall"))
  # Check thresh.method
  thresh.method <- match.arg(thresh.method,c("absolute","density","p-value"))
  # Check adjust
  if(thresh.method!="p-value" & !missing(adjust)) stop("adjust can be set only for thresh.method equal to p-value")
  if(thresh.method=="pvalue"){
    adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  }
  # Check thresh
  if(!is.numeric(thresh) | thresh<0 | thresh>1) stop("thresh must number in range [0,1]")
  
  # Add 1 if are present zeros
  ifelse(any(mg@data==0), x<-mg@data+1, x<-mg@data)
  
  if(thresh.method=="absolute"){
    adj <- cor(mgnet::clr(x),method=cor.method)
    adj <- adj * (adj>=thresh)
  } else if(thresh.method=="density"){
    adj <- mgnet::adjacency_edge_density(x=mgnet::clr(x),method=cor.method,th=thresh)
  } else if(thresh.method=="p-value"){
    adj <- mgnet::adjacency_p_adjust(x=mgnet::clr(x),method=cor.method,adjust=adjust,alpha=thresh)
  }
  
  return(mgnet(mg=mg,adj=adj))
})
#' @rdname make_mgnet-methods
#' @aliases make_mgnet,mg,character,character,missing,numeric
setMethod("make_mgnet", c("mg","character","character","missing","numeric"),
          function(mg, cor.method, thresh.method, adjust, thresh){
            
            # Check mg
            if(!isa(mg,"mg")) stop("mg must belong to mg class")
            # Check cor.method
            cor.method <- match.arg(cor.method,c("pearson","spearman","kendall"))
            # Check thresh.method
            thresh.method <- match.arg(thresh.method,c("absolute","density"))
            # Check thresh
            if(!is.numeric(thresh) | thresh<0 | thresh>1) stop("thresh must number in range [0,1]")
            
            # Add 1 if are present zeros
            ifelse(any(mg@data==0), x<-mg@data+1, x<-mg@data)
            
            if(thresh.method=="absolute"){
              adj <- cor(mgnet::clr(x),method=cor.method)
              adj <- adj * (adj>=thresh)
            } else if(thresh.method=="density"){
              adj <- mgnet::adjacency_edge_density(x=mgnet::clr(x),method=cor.method,th=thresh)
            }
            
            return(mgnet(mg=mg,adj=adj))
          })
#' @rdname make_mgnet-methods
#' @aliases make_mgnet,list,character,character,character,numeric
setMethod("make_mgnet", c("list","character","character","character","numeric"),
          function(mg, cor.method, thresh.method, adjust, thresh){
            lapply(mg, 
                   selectMethod(f="make_mgnet",signature=c("mg","character","character","character","numeric")),
                   cor.method=cor.method, thresh.method=thresh.method,
                   adjust=adjust, thresh=thresh)
          })
#' @rdname make_mgnet-methods
#' @aliases make_mgnet,list,character,character,missing,numeric
setMethod("make_mgnet", c("list","character","character","missing","numeric"),
          function(mg, cor.method, thresh.method, adjust, thresh){
            lapply(mg, 
                   selectMethod(f="make_mgnet",signature=c("mg","character","character","missing","numeric")),
                   cor.method=cor.method, thresh.method=thresh.method,
                   thresh=thresh)
          })
>>>>>>> main
