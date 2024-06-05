# CLR !!!
#------------------------------------------------------------------------------#
#' Centered Log-Ratio Transformation
#'
#' @description
#' Performs the Centered Log-Ratio (CLR) transformation on a given dataset (or at least a single vector). 
#' This transformation is a standard approach for analyzing compositional data, 
#' which consists of proportions or percentages summing up to a constant sum. 
#' It transforms the data into a space where standard statistical methods can be applied. 
#' Rows represent samples, and columns represent components or taxa.
#'
#' @param X A numeric matrix with all elements greater than 0, where rows are samples 
#' and columns are components or taxa.
#'
#' @return A numeric matrix of the same dimension as `X`, representing the clr-transformed data.
#'
#' @export
clr <- function(X) {
  
  if( !(is.matrix(X) | is.vector(X)) ) stop("X must be a matrix or a vector")
  if(!is.numeric(X) | any(X<=0)) stop("X must be numeric with all elements greater than 0")

  if(is.null(dim(X))){
    return(log(X) - mean(log(X)))
  } else {
    ref <- apply(X, 1, function(x) mean(log(x)) )
    return(as.matrix(log(X) - ref))
  }
  
}


# ICLR !!!
#------------------------------------------------------------------------------#
#' Interquartile Centered Log-Ratio Transformation
#'
#' @description
#' Performs the Interquartile Centered Log-Ratio (iclr) transformation on a given dataset. 
#' This transformation is designed to reduce the influence of highly variable components 
#' by normalizing the data using a subset of components with lower variability.
#' The transformation is applied to compositional data where rows represent samples 
#' and columns represent components or taxa.
#'
#' @param X A numeric matrix with all elements greater than 0, where rows are samples 
#' and columns are components or taxa.
#'
#' @return A numeric matrix of the same dimension as `X`, representing the iclr-transformed data.
#'
#' @export
#' @importFrom stats var quantile
iclr <- function(X) {
  if (!is.numeric(X) || any(X <= 0)) {
    stop("X must be numeric with all elements greater than 0.")
  }
  
  # Ensure X is a matrix, necessary for operations below
  if (is.vector(X)) {
    X <- matrix(X, nrow = 1)
  }
  
  # CLR transformation
  X_log <- log(X)
  geomean <- rowMeans(X_log, na.rm = TRUE)
  X_clr <- sweep(X_log, 1, geomean, FUN = "-")
  
  # Calculate the variance of CLR-transformed values and find the interquartile range
  clr_var <- apply(X_clr, 2, var)
  iqr_limits <- quantile(clr_var, probs = c(0.25, 0.75), na.rm = TRUE)
  
  # Determine the invariant set (features within the interquartile range of variability)
  invariant_set <- clr_var >= iqr_limits[1] & clr_var <= iqr_limits[2]
  
  # Refine the CLR transformation using only the invariant set
  if (any(invariant_set)) {
    geomean_inv <- rowMeans(X_log[, invariant_set], na.rm = TRUE)
    X_iqclr <- sweep(X_log, 1, geomean_inv, FUN = "-")
  } else {
    # If no invariant set is found, return the original CLR transformation
    X_iqclr <- X_clr
  }
  
  return(X_iqclr)
}