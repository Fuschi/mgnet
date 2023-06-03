################################################################################
################################################################################
# ZERO_DEALING
################################################################################
################################################################################
#' Simplest way to replace zero.
#' 
#' @description This function implement the simplest replacement of the zeros 
#' in a compositional data sets.
#'
#' @param X compositional matrix.
#' @param mar integer indicating the margin of samples (rows=1 or cols=2).
#' @param type character indicating the zero replacement strategy with 
#' possible values "add" or "sub". If type is equal to "add" a pseudocount is 
#' added to all elements of X instead if it is set "sub" the  only the zeros are
#' replaced.
#' 
#' @export
zero_dealing <- function(X, mar=1, type="plus"){
  
  if(!is.numeric(X) | !is.matrix(X)) stop("X must be numeric matrix")
  if(any(X<0)) stop("Find negative elements in X")
  if(!(mar%in%c(1,2))) stop("mar has only 1 or 2 as possibles values")
  type <- match.arg(type,c("sum","sub"))
  
  if(all(round(X)-X==0)){
    dl <- matrix(.65,nrow=nrow(X),ncol=ncol(X))
  } else {
    dl <- apply(X, mar, function(x)min(x[x>0]))
    ifelse(mar==1,
           dl <- replicate(ncol(X),dl),
           dl <- t(replicate(nrow(X),dl)))
    dl <- .65*dl
  }
  
  if(type=="plus"){
    Y <- X + dl
  } else {
    Y <- X + dl*(X==0)
  }
  
  return(Y)
}
################################################################################
################################################################################
# END ZERO_DEALING REPLACEMENT
################################################################################
################################################################################