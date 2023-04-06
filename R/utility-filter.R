#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Prevalence
#' 
#' @description Return the portion of the elements greater than zero of a vector
#' or a matrix. 
#' 
#' @param X vector or matrix with all elements >=0.
#' @param mar a vector giving the subscripts which the function will be applied over; 1 indicates rows, 2 indicates columns, 0 for all elements.
#' 
#' @export
prevalence <- function(X, mar=2){
  if(!(is.vector(X) | is.matrix(X))) stop("X must be a matrix or a vector")
  if(!is.numeric(X)) stop("X must be numeric")
  if(any(X<0)) stop("All elements of X must be >=0")
  if(!is.numeric(mar)) stop("mar must be numeric")
  if(!any(mar%in%c(0,1,2))) stop("mar must be equal at one of the following choices 1,2 or 0")
  
  if(is.vector(X)){
    return(sum(X>0)/length(X))
  } else if (mar==1 | mar==2){
    return(apply(X,mar,function(x) sum(x>0)/length(x), simplify=T))
  } else {
    return(sum(X>0)/length(X))
  }
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Quantile Non Zero Elements 
#' 
#' @description Return the choosen quantile of the elements greater than zero from a vector
#' or a matrix. 
#' 
#' @param X vector or matrix with all elements >=0.
#' @param mar a vector giving the subscripts which the function will be applied over; 1 indicates rows, 2 indicates columns, 0 for all elements.
#' @param probs numeric vector of probabilities with values in \[0,1\] with default
#' set to 0.5, thus the median values (see \code{\link{quantile}}).
#' @param ... additional arguments of quantile.
#' 
#' @importFrom stats quantile
#' @export
quantile_nozero <- function(X, probs=.5, mar=2, ...){
  if(!(is.vector(X) | is.matrix(X))) stop("X must be a matrix or a vector")
  if(!is.numeric(X)) stop("X must be numeric")
  if(any(X<0)) stop("All elements of X must be >=0")
  if(!is.numeric(mar)) stop("mar must be numeric")
  if(!any(mar%in%c(0,1,2))) stop("mar must be equal at one of the following choices 1,2 or 0")
  if(!is.numeric(probs)) stop("probs must be numeric")
  if(any( probs<0 | probs>1 )) stop("the elements of probs must be in range [0,1]")
  
  if(is.vector(X)){
    
    if(any(X>0)){
      return(quantile(X[X>0],probs=probs,...))
    } else {
      return(0)
    }
    
  } else if (mar==1 | mar==2){
    
    return(apply(X,mar,function(x) ifelse(any(x>0),
                                          quantile(x[x>0],probs=probs,...),
                                          0)))
    
  } else {
    
    if(any(c(X)>0)){
      return(quantile(c(X)[c(X)>0],probs=probs,...))
    } else {
      return(0)
    }
    
  }
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Shannon alpha diversity
#' 
#' @description Returns the Shannon alpha diversity from the margins of a matrix.
#'  
#' 
#' @param X vector or matrix with all elements >=0.
#' @param mar a vector giving the subscripts which the function will be applied over; 1 indicates rows, 2 indicates columns, 0 for all elements.
#' @param norm logical elements indicates if the entropy must be normalize to 1 dividing it
#' by the logarithm of the dimension (default TRUE).
#' 
#' @export
shannon_diversity <- function(X, mar=2, norm=TRUE){
  if(!(is.vector(X) | is.matrix(X))) stop("X must be a matrix or a vector")
  if(!is.numeric(X)) stop("X must be numeric")
  if(any(X<0)) stop("All elements of X must be >=0")
  if(!is.numeric(mar)) stop("mar must be numeric")
  if(!any(mar%in%c(0,1,2))) stop("mar must be equal at one of the following choices 1,2 or 0")
  if(!is.logical(norm)) stop("norm must be logical")
  
  if(is.vector(X)){
    p <- X/sum(X)
    e <- -sum(p*log(p), na.rm=T)
    if(norm) e <- e/log(length(X))
  } else {
    
    if(mar==1){
      p <- X/rowSums(X)
      e <- -rowSums(p*log(p), na.rm=T)
      if(norm) e <- e/log(ncol(X))
    } else if(mar==2){
      p <- X/colSums(X)
      e <- -colSums(p*log(p), na.rm=T)
      if(norm) e <- e/log(nrow(X))
    } else {
      p <- c(X)/sum(c(X))
      e <- -sum(p*log(p), na.rm=T)
      if(norm) e <- e/log(length(X))
    } 
    
  }
  
  return(e)
}