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
#' @param type character indicating the zero replacement strategy. Only two
#' possible values are available, plus and subs.
#' 
#' @export
zero_dealing <- function(X, mar=1, type="plus"){
  
  if(!is.numeric(X) | !is.matrix(X)) stop("X must be numeric matrix")
  if(any(X<0)) stop("Find negative elements in X")
  if(all(X>0)) warning("There are no zeros in X")
  if(!(mar%in%c(1,2))) stop("mar has only 1 or 2 as possibles values")
  type <- match.arg(type,c("plus","subs"))
  
  if(all(round(X)-X==0)){
    dl <- matrix(1,nrow=nrow(X),ncol=ncol(X))
  } else {
    dl <- apply(X, mar, function(x)min(x[x>0]))
    ifelse(mar==1,
           dl <- replicate(ncol(X),dl),
           dl <- t(replicate(nrow(X),dl)))
  }
  
  if(type=="plus"){
    Y <- X + .65*dl
  } else {
    Y <- X + .65*dl*(X==0)
  }
  
  return(Y)
}
################################################################################
################################################################################
# END ZERO_DEALING REPLACEMENT
################################################################################
################################################################################