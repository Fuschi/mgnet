#' Get Adjacency with Corrected P-values
#'
#' @description
#' Generates an adjacency matrix from abundance data or a pre-computed correlation matrix,
#' applying a threshold on the correlation matrix based on adjusted p-values with multiple
#' test corrections. This function is useful for constructing networks from correlation matrices
#' by considering statistically significant correlations.
#'
#' @param x A numeric matrix of data, where each row represents a sample and each column represents a variable/taxon.
#'          This parameter is required if `r` and `n` are not provided.
#' @param r A pre-computed correlation matrix. If provided, `n` must also be specified.
#' @param n The number of observations used to compute the correlation matrix. Required if `r` is provided.
#' @param method Character string specifying the correlation method ("pearson", "spearman", "kendall").
#'               Default is "pearson".
#' @param adjust Character string specifying the method for p-value adjustment ("holm", "hochberg", "hommel",
#'               "bonferroni", "BH", "BY", "fdr", "none"). Default is "holm".
#' @param alpha Numeric value specifying the significance level for the adjusted p-values. Default is 0.05.
#'
#' @return An adjacency matrix derived from the correlation matrix, where elements are set to zero
#'         if their corresponding adjusted p-value is greater than `alpha`.
#'
#' @export
#' @importFrom stats p.adjust cor
adjacency_p_adjust <- function(x = NULL, r = NULL, n = NULL, method = "pearson", adjust = "holm", alpha = 0.05) {
  
  # CHECKS
  #----------------------------------------------------------------------------#

  if (!is.null(x) && (!is.null(r) || !is.null(n))) {
    stop("Provide either 'x' OR both 'r' and 'n', not a combination.")
  }
  
  if (!is.null(r) && is.null(n)) {
    stop("'r' and 'n' must be provided together.")
  }
  
  if(!is.null(x) && (!is.matrix(x) || !is.numeric(x))) {
    stop("'x' must be a numeric matrix.")
  }
  
  if(!is.null(r)) {
    if(!is.matrix(r) || !is.numeric(r) || !isSymmetric(r) || any(abs(r) > 1)) {
      stop("'r' must be a numeric symmetric matrix with elements in range [-1,1].")
    }
    if(!is.numeric(n) || n <= 0 || round(n) != n) {
      stop("'n' must be a positive integer.")
    }
  }
  
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  adjust <- match.arg(adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # Calculate correlation matrix if not provided
  if (is.null(r)) {
    r <- cor(x, use = "pairwise.complete.obs", method = method)
    n <- nrow(x) 
  }
  
  # Calculate p-values from the correlation matrix
  p_vals <- matrix(nrow = ncol(r), ncol = ncol(r), dimnames = list(colnames(r), colnames(r)))
  for (i in 1:(ncol(r) - 1)) {
    for (j in (i + 1):ncol(r)) {
      test <- cor.test(r[, i], r[, j], method = method)
      p_vals[i, j] <- test$p.value
      p_vals[j, i] <- test$p.value
    }
  }
  
  # Adjust p-values
  if (adjust != "none") {
    p_vals[upper.tri(p_vals, diag = FALSE)] <- p.adjust(p_vals[upper.tri(p_vals, diag = FALSE)], method = adjust)
    p_vals[lower.tri(p_vals, diag = FALSE)] <- t(p_vals)[lower.tri(p_vals, diag = FALSE)]
  }
  
  # Create adjacency matrix
  adjacency_matrix <- (r * (p_vals <= alpha))
  diag(adjacency_matrix) <- 0  
  
  return(adjacency_matrix)
}


#' Generate Adjacency Matrix by Edge Density Threshold
#'
#' @description
#' Creates an adjacency matrix from abundance data or a pre-computed correlation matrix
#' by applying a threshold based on the desired edge density of the resulting network.
#'
#' @param x A numeric matrix of data, where each row represents a sample and each column represents a variable.
#'          Used for computing the correlation matrix if `r` is not provided. Either `x` or `r` must be provided.
#' @param r A pre-computed correlation matrix. If provided, `x` should not be specified.
#' @param method Character string specifying the correlation method ("pearson", "spearman", "kendall").
#'               Used only if `x` is provided. Default is "pearson".
#' @param th Numeric value specifying the desired edge density threshold in range \[0,1\]. Default is 0.05.
#'
#' @return An adjacency matrix where only correlations that meet the edge density threshold are retained.
#' 
#' @importFrom stats cor
#' @export
adjacency_edge_density <- function(x=NULL, r=NULL, method="pearson", th=.05){
  
  if(!is.null(x) && !is.null(r)) {
    stop("Provide only 'x' OR 'r', not both.")
  }
  
  if(!is.null(x) && (!is.matrix(x) || !is.numeric(x))) {
    stop("'x' must be a numeric matrix.")
  }
  
  if(!is.null(r)) {
    if(!is.matrix(r) || !is.numeric(r) || !isSymmetric(r) || any(abs(r) > 1)) {
      stop("'r' must be a numeric symmetric matrix with elements in range [-1,1].")
    }
  }
  
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  
  if(!is.numeric(th) || th < 0 || th > 1) {
    stop("'th' must be a number in range [0,1].")
  }
  
  # Calculate correlation matrix if 'x' is provided
  if(!is.null(x)) {
    r <- cor(x, method=method)
  }
  
  diag(r) <- 0
  
  # Take correlation upper triangular values as vector
  rtriu <- r[upper.tri(r)]
  
  # Sort Absolute values
  rtriu.sort <- sort(abs(rtriu),decreasing=TRUE)
  
  # Take the Absolute Min Values
  th.position <- round(th*(ncol(r)*(ncol(r)-1))*.5)
  th.value <- rtriu.sort[th.position]
  
  adj <- r * (abs(r)>=th.value)
  return(adj)
}


#' Generate Adjacency Matrix by Absolute Value Threshold
#'
#' This function creates an adjacency matrix by applying a threshold to the absolute 
#' values of the correlations in a correlation matrix. Elements of the correlation 
#' matrix greater than or equal to the threshold are preserved, while others are set to zero.
#' This approach is useful in network analysis, where an edge between two nodes (variables)
#' is considered significant if their correlation's absolute value exceeds the specified threshold.
#'
#' @param x Optional; a numeric matrix of data from which to calculate the correlation matrix. 
#' Provide either 'x' or 'r', but not both. If 'x' is provided, 'r' should be NULL.
#' @param r Optional; a pre-computed numeric symmetric correlation matrix with elements in 
#' the range \[-1,1\]. Provide either 'x' or 'r', but not both. If 'r' is provided, 'x' should be NULL.
#' @param method Character; the method for computing correlation if 'x' is provided. 
#' Acceptable values are "pearson" (default), "spearman", or "kendall".
#' @param th Numeric; the threshold for determining significance in the adjacency matrix. 
#' Must be a number in the range \[0,1\]. Correlations with absolute values greater than or equal 
#' to this threshold are considered significant.
#'
#' @return A numeric adjacency matrix with the same dimensions as the correlation matrix 'r'.
#' Elements are set to the value of the correlation if the absolute value of the correlation 
#' is greater than or equal to 'th', and set to zero otherwise.
#' 
#' @importFrom stats cor
#' @export
adjacency_absolute_value <- function(x=NULL, r=NULL, method="pearson", th=.3){
  
  if(!is.null(x) && !is.null(r)) {
    stop("Provide only 'x' OR 'r', not both.")
  }
  
  if(!is.null(x) && (!is.matrix(x) || !is.numeric(x))) {
    stop("'x' must be a numeric matrix.")
  }
  
  if(!is.null(r)) {
    if(!is.matrix(r) || !is.numeric(r) || !isSymmetric(r) || any(abs(r) > 1)) {
      stop("'r' must be a numeric symmetric matrix with elements in range [-1,1].")
    }
  }
  
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  
  if(!is.numeric(th) || th < 0 || th > 1) {
    stop("'th' must be a number in range [0,1].")
  }
  
  # Calculate correlation matrix if 'x' is provided
  if(!is.null(x)) {
    r <- cor(x, method=method)
  }
  
  diag(r) <- 0
  adj <- r * (abs(r)<=th)
  return(adj)
}