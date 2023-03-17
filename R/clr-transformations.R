################################################################################
################################################################################
# CLR
################################################################################
################################################################################
#' Centered Log-Ratio Transformation
#' 
#' @description It calculates the centered log-ratio transformation of X.
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
# INTERQUARTILE CLR
################################################################################
################################################################################
#' Robust Centered Log-Ratio Transformation
#'
#' @description It calculates the interquantile centered log-ratio transformation of X.
#'
#' @param X (Required) numeric matrix or vector with all elements greater than or equal to 0.
#' @param mar (Optional) Integer giving the dimension where the function will be applied;
#' 1 for rows and 2 for columns (default 1).
#' 
#' @importFrom stats var quantile
#' 
#' @export
iqclr <- function(X, mar=1){
  
  if( !(is.matrix(X)) ) stop("X must be a matrix")
  if(!is.numeric(X) | any(X<=0)) stop("X must be numeric with all elements greater than 0")
  if(!(mar%in%c(1,2))) stop("mar has as possible values only 1 and 2.")
  
  Xclr <- mgnet::clr(X, mar=mar)
  
  # find iqlr denom
  clr_var=apply(Xclr, c(1,2)[mar], stats::var)
  qts <- stats::quantile(clr_var, na.rm=T, probs=seq(0,1,by=0.25))
  invariant.set <- which( clr_var>=qts[2] & clr_var<=qts[4] )
  
  if(mar==1){
    Xinv <- X[,invariant.set]
  } else {
    Xinv <- X[invariant.set,]
  }
  
  ref <- apply(Xinv, mar, function(x) mean(log(x)) )
  return(as.matrix(log(X) - ref))
  
}
################################################################################
################################################################################
# END CLR
################################################################################
################################################################################