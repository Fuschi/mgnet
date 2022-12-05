################################################################################
################################################################################
# CLR
################################################################################
################################################################################
#' Centered Log-Ratio Transformation
#'
#' @param X (Required) numeric matrix or vector with all elements greater than or equal to 0.
#' @param mar (Optional) Integer giving the dimension where the function will be applied;
#' 1 for rows and 2 for columns (default 1).
#' 
#' @export
clr <- function(X, mar=1){
  
  if( !(is.matrix(X) | is.vector(X)) ) stop("X must be a matrix or a vector")
  if(!is.numeric(X) | any(X<=0)) stop("X must be numeric with all elements greater than 0")
  if(!(mar%in%c(1,2))) stop("mar has as possible values only 1 and 2.")
  
  if(is.null(dim(X))){
    return(log(X) - mean(log(X)))
  } else {
    ref <- apply(X, mar, function(x) mean(log(x)) )
    return(as.matrix(log(X) - ref))
  }
}
################################################################################
################################################################################
# END CLR
################################################################################
################################################################################




################################################################################
################################################################################
# ADJACENCY P-ADJUST
################################################################################
################################################################################
#' Get adjacency with corrected p-values
#' 
#' @description User wrapper of \code{\link{corr.test}} or \code{\link{corr.p}}
#' to get the adjacency from abundances or correlation matrix. The function apply a
#' threshold on correlation matrix using the adjusted p-value with multiple test
#' corrections.
#' 
#' @param x matrix of data.
#' @param r correlation matrix.
#' @param n number of observation.
#' @param method (Optional) correlation type, possible choices are "pearson",
#' "spearman" or "kendall" (default spearman).
#' @param adjust (Optional) multiple test choosen, possible choices are "holm",
#' "hochberg","hommel","bonferroni","BH","BY","fdr","none".
#' @param alpha (Optional) level of confidence interval.
#' @param keep.psych (Optional) logical paramater which indicates if return the
#' results of psych results (default FALSE).
#' 
#' 
#' @importFrom psych corr.test corr.p
#' @export
adjacency_p_adjust <- function(x=NULL,r=NULL,n=NULL,
                           method="pearson",adjust="holm",alpha=.05,
                           keep.psych=FALSE){
  
  # Check the logical assignment of the parameters.
  if(!is.null(x) & (!is.null(r) | !is.null(n))) stop("x, r and n can't be not null togheter. The function permits the cases the abundances x OR the correlation matrix r with the number of observation n.")
  if((!is.null(r) & is.null(n)) | (is.null(r) & !is.null(n)) ) stop("To compute the adjusted p from the correlation matrix r then the observation number n must also be assigned (and vice versa)")
  
  # Check parameters properties
  method <- match.arg(method, c("pearson","spearman","kendall"))
  adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  if(!is.numeric(alpha)) stop("alpha must be a number")
  if(alpha<0 | alpha>1) stop("alpha must be a number in range [0,1]")
  if(!is.logical(keep.psych)) stop("keep.psych must be logical")
  if(!is.null(x)){
    if(!is.numeric(x) | !is.matrix(x)) stop("x must be a numeric matrix")
  }
  if(!is.null(r)){
    if(!is.numeric(r) | !is.matrix(r) | !isSymmetric(r)) stop("x must be a numeric symmetric matrix")
    if(any(abs(r)>1)) stop("elements of r must be in range [-1,1]")
    if(n<=0 | round(n)!=n) stop("n must be positive integer number")
  }
  
  
  # Get adjusted p-values
  if(!is.null(x)){
    res.psych <- psych::corr.test(x=x,use="pairwise",method=method,adjust=adjust,
                                  alpha=alpha,ci=FALSE)
    padj <- res.psych$p;
    padj[upper.tri(padj)] <- res.psych$p.adj
    padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]
    diag(padj) <- 1
    adj <- res.psych$r*(padj<=alpha)
    
  } else {
    
    res.psych <- psych::corr.p(r=r,n=n,adjust=adjust,
                                  alpha=alpha,ci=FALSE)
    padj <- res.psych$p;
    padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]
    diag(padj) <- 1
    adj <- res.psych$r*(padj<=alpha)

  }
  
  # return
  if(keep.psych){return(list("adj"=adj,
                "psyc"=res.psych))
  } else {return(adj)}
  
}
################################################################################
################################################################################
# END ADJACENCY P-ADJUST
################################################################################
################################################################################




################################################################################
################################################################################
# ADJACENCY EDGE-DENSITY
################################################################################
################################################################################
#' Get adjacency making an edge density threshold
#' 
#' @description Retrieves the adjacency matrix starting from the abundance 
#' or the correlation matrix. The function uses the \code{\link{cor}} 
#' function from stats.
#' 
#' @param x matrix of data.
#' @param r correlation matrix.
#' @param method (default "pearson") correlation type, possible choices are "pearson",
#' "spearman" or "kendall" (default spearman).
#' @param th (default 0.05) edge density choosen.
#' 
#' @importFrom stats cor
#' 
#' @export
adjacency_edge_density <- function(x=NULL,r=NULL, method="pearson",th=.05){
  
  # Check the logical assignment of the parameters.
  if(!is.null(x) & !is.null(r)) stop("x, r and n can't be not null togheter. The function permits the cases of abundances x OR the correlation matrix r.")
  
  # Check parameters properties
  method <- match.arg(method, c("pearson","spearman","kendall"))
  if(!is.numeric(th) | th<0 | th>1) stop("th must be a number in range [0,1]")
  if(!is.null(x)){
    if(!is.numeric(x) | !is.matrix(x)) stop("x must be a numeric matrix")
  }
  if(!is.null(r)){
    if(!is.numeric(r) | !is.matrix(r) | !isSymmetric(r)) stop("x must be a numeric symmetric matrix")
    if(any(abs(r)>1)) stop("correlation elements of r must be in range [-1,1]")
  }
  
  # Get Correlation from x
  if(!is.null(x)) r <- stats::cor(x,method=method)
  
  diag(r) <- 0
  
  # Take correlation upper triangular values as vector
  rtriu <- r[upper.tri(r)]
  
  # Sort Absolute values
  rtriu.sort <- sort(abs(rtriu),decreasing=TRUE)
  
  # Take the Absolute Min Values
  th.value <- rtriu.sort[round(length(rtriu.sort)*th)]
  
  adj <- r * (abs(r)>=th.value)
  return(adj)
  
}
################################################################################
################################################################################
# END ADJACENCY EDGE-DENSITY
################################################################################
################################################################################
