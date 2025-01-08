################################################################################
################################################################################
# ZERO_DEALING
################################################################################
################################################################################
#' Zero Replacement in Compositional Data
#'
#' @description
#' This function implements the simplest strategies for zero replacement in compositional datasets. 
#' Zero replacement is a critical step in preprocessing compositional data, as it 
#' allows for the application of compositional data analysis techniques which assume 
#' non-zero data (like for log-ratio transformetions). 
#'
#' The "unif" method replaces zeros with random values drawn from a uniform distribution, 
#' with bounds set to be a small fraction of the minimum non-zero value within each sample. 
#'
#' The "const" method, on the other hand, replaces zeros with a constant value, 
#' specifically 65% of the minimum non-zero value found in each sample. 
#'
#' @section Discussion: 
#' It's important to note that these methods are among the simplest and most 
#' straightforward approaches to zero replacement. While they provide a practical
#' solution for dealing with zeros in compositional data, especially in datasets 
#' with a large number of zeros, they may not always be the best choice for all 
#' types of compositional analyses. 
#'
#' @param X A numeric matrix representing the compositional data, where each row represents a sample and each column represents a variable/taxon.
#' @param method A character string indicating the zero replacement strategy: "unif" for uniform random values or "const" for a constant value replacement.
#'
#' @return A numeric matrix with the same dimensions as `X`, where zeros have been replaced according to the specified method.
#' @export
#' @references
#' Lubbe S. et al, 2021. Comparison of zero replacement strategies for compositional data with large numbers of zeros. Chemometrics and Intelligent Laboratory Systems. 10.1016/j.chemolab.2021.104248
#'
#' @seealso
#' The following packages provide a range of methodologies for addressing zeros 
#' in compositional data, from statistical imputation techniques to advanced machine 
#' learning models, offering suitable options for different types of analysis and data complexity.
#' \itemize{
#'  \item zCompositions-package
#'  \item robCompositions-package
#'  \item deepImp-package
#' }
#' 
zero_dealing <- function(X, method = "const"){
  
  # Checks
  if(!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if(any(X < 0)) stop("X contains negative elements. Only non-negative values are allowed.")
  method <- match.arg(method, c("unif", "const"))
  if( any(rowSums(X)==0) ) stop("A sample with all elements equal to 0 was found.")
  
  # Check if X contains only integer values
  is_integer_matrix <- all(X == floor(X))
  
  if(is_integer_matrix & method == "unif"){
    
    idx_zero <- which(X==0)
    X[idx_zero] <- stats::runif(n=length(idx_zero), min = .065, max = .65)
    return(X)
    
  } else if (is_integer_matrix & method == "const") {
    
    X[X==0] <- .65
    return(X)
    
  } else if (!is_integer_matrix & method == "unif") {
    
    # zero uniform replacement on a single vector/sample
    zero_unif <- function(x){
      minx <- min(x[x>0])
      idx_zero <- which(x==0)
      x[idx_zero] <- stats::runif(length(idx_zero), min = .1*minx, max = minx)
      return(x)
    }
    
    X <- apply(X,1,zero_unif)
    return(t(X))
    
  } else if (!is_integer_matrix & method == "const") {
    
    # zero uniform replacement on a single vector/sample
    zero_const <- function(x){
      minx <- min(x[x>0])
      idx_zero <- which(x==0)
      x[idx_zero] <- .65*minx
      return(x)
    }
    
    X <- apply(X,1,zero_const)
    return(t(X)) 
    
  } else {
    
    stop("I'm lost... Why are we still here? Just to suffer? cit. 'Kaz'\nJoke aside, you shouldn't be here")
    
  }
  
}
################################################################################
################################################################################
# END ZERO_DEALING REPLACEMENT
################################################################################
################################################################################


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


#' CLR Transformation with Zero Replacement
#'
#' This function performs a centered log-ratio (CLR) or isometric CLR (iCLR) transformation 
#' on an abundance matrix. It includes options for zero replacement strategies, allowing
#' for the substitution of zero values with either a constant or a uniform random value.
#'
#' @param X A numeric matrix representing abundance or relative abundance data. The matrix should not be empty.
#' @param clr_variant A character string specifying the variant of log-ratio transformation to apply.
#'        Either `"clr"` (centered log-ratio) or `"iclr"` (isometric centered log-ratio).
#'        Default is `"clr"`.
#' @param zero_strategy A character string specifying the strategy for replacing zero values in the matrix.
#'        Either `"const"` to replace zeros with a constant value or `"unif"` to replace zeros with a uniform random value.
#'        Default is `"const"`.
#'
#' @return A matrix that has been transformed using the CLR or iCLR method with zero handling applied.
#'
#' @details
#' The function first replaces zero values in the matrix according to the chosen `zero_strategy`. 
#' It then applies either the CLR or iCLR transformation. 
#' The CLR transformation transforms the data into log-ratio form by dividing each value 
#' by the geometric mean of its corresponding row, followed by taking the log. 
#' The iCLR transformation is a variation of CLR that maintains isometry in the transformed space.
#' 
#' @references
#' Lubbe S. et al, 2021. Comparison of zero replacement strategies for compositional data with large numbers of zeros. Chemometrics and Intelligent Laboratory Systems. 10.1016/j.chemolab.2021.104248
#'
#' @seealso 
#' \code{\link[=clr]{clr}} for standard centered log-ratio transformation,
#' \code{\link[=iclr]{iclr}} for isometric centered log-ratio transformation,
#' and \code{\link[=zero_dealing]{zero_dealing}} for zero handling strategies.
#'
#'
#' @export
clr_zero_handle <- function(X, clr_variant = c("clr", "iclr"), zero_strategy = c("const", "unif")) {
  # Validate arguments
  clr_variant <- match.arg(clr_variant, c("clr", "iclr"))
  zero_strategy <- match.arg(zero_strategy, c("const", "unif"))
  
  # Handle zero values in the matrix based on the selected strategy
  X_nozero <- zero_dealing(X = X, method = zero_strategy)
  
  # Apply clr or iclr transformation based on the clr_variant argument
  if (clr_variant == "clr") {
    return(clr(X_nozero))
  } else {
    return(iclr(X_nozero))
  }
}
