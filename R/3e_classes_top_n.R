# TOP_N_SAMPLE
#------------------------------------------------------------------------------#
#' Extract Top N Samples Based on Aggregate Metrics
#'
#' Selects the top `n` samples from an `mgnet` object based on an aggregate metric calculated from 
#' specified data fields such as abundance, relative abundance, or normalized abundance, with
#' the chosen metric, which can be customized.
#'
#' @param object An `mgnet` object.
#' @param n Integer specifying the number of top samples to retrieve based on the calculated metric.
#' @param field Character string specifying the data field to use for computing the metric.
#'        Acceptable values are `"abundance"`, `"rel_abundance"`, or `"norm_abundance"`.
#'        The field must not be empty in the `mgnet` object.
#' @param rank Optional character string specifying the taxonomic rank to consider when calculating metrics.
#'        If provided, the computation aggregates data at this taxonomic rank within the `lineage` data.
#'        Note: Aggregating data at different taxonomic levels may sum taxa classified together,
#'        potentially affecting the interpretability of normalized matrices.
#' @param order_fun Function used to compute the aggregate metric across the selected field. The default 
#'        is `sum`, but other functions like `mean` or `median` can be used to tailor the analysis.
#' @param decreasing Logical indicating if sorting should be in decreasing order. When `TRUE` (default),
#'        samples with the highest metric values are considered top.
#'
#' @return For `mgnet` objects, a character vector of sample IDs representing the top `n` samples.
#'         For `mgnetList` objects, returns a list where each element is a character vector of top `n`
#'         samples for each respective `mgnet`.
#' @note Ensure that the specified field is not empty in the `mgnet` object before calling this function. 
#'       Failing to do so will result in an error. Aggregating data at different taxonomic levels
#'       is the result of summing up the same taxa together and not always it has sense for norm_abundance.
#' @export
#' @name top_n_samples
#' @aliases top_n_samples,mgnet-method top_n_samples,mgnetList-method
setGeneric("top_n_samples", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {standardGeneric("top_n_samples")})

setMethod("top_n_samples", "mgnet", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {

  if (!field %in% c("abundance", "rel_abundance", "norm_abundance")) {
    stop("Invalid field specified; choose from 'abundance', 'rel_abundance', 'norm_abundance'.")
  } else if(field == "abundance" & length(object@abundance) == 0) {
    stop("abundance slot cannot be empty if field is setted to 'abundance'")
  } else if(field == "rel_abundance" & length(object@rel_abundance) == 0) {
    stop("rel_bundance slot cannot be empty if field is setted to 'rel_abundance'")
  } else if(field == "norm_abundance" & length(object@norm_abundance) == 0) {
    stop("norm_abundance slot cannot be empty if field is setted to 'norm_abundance'")
  }
  
  if (!is.null(rank) && !rank %in% colnames(object@lineage)) {
    stop("Specified rank not found in lineage data.")
  } 
  
  # Extract the data matrix
  if(is.null(rank)){
    
    if(field == "abundance"){
      data_matrix <- abundance(object)
    } else if(field == "rel_abundance"){
      data_matrix <- rel_abundance(object)
    } else if(field == "norm_abundance"){
      data_matrix <- norm_abundance(object)
    }
    
  } else {
    
    if(field == "abundance"){
      data_matrix <- abundance(object, rank = rank)
    } else if(field == "rel_abundance"){
      data_matrix <- rel_abundance(object, rank = rank)
    } else if(field == "norm_abundance"){
      data_matrix <- norm_abundance(object, rank = rank)
    }
    
  }
  
  
  # Apply the order function to sum up or calculate the desired metric across columns (taxa)
  sample_scores <- apply(data_matrix, 1, order_fun)
  
  # Order by the computed scores
  if (decreasing) {
    sample_indices <- order(sample_scores, decreasing = TRUE)
  } else {
    sample_indices <- order(sample_scores, decreasing = FALSE)
  }
  
  # Subset to top n
  top_n_indices <- sample_indices[1:n]
  
  return(sample_id(object)[top_n_indices])
  
})

setMethod("top_n_samples", "mgnetList", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {
  
  result <- sapply(object, FUN = \(x){
    top_n_samples(x, n = n, field = field, rank = rank, 
                  order_fun = order_fun, decreasing = decreasing)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(result)
})


# TOP_N_SAMPLE
#------------------------------------------------------------------------------#
#' Extract Top N Taxa Based on Aggregate Metrics
#'
#' Selects the top `n` taxa from an `mgnet` object based on an aggregate metric calculated from 
#' specified data fields such as abundance, relative abundance, or normalized abundance, with
#' the chosen metric, which can be customized.
#'
#' @param object An `mgnet` object.
#' @param n Integer specifying the number of top taxa to retrieve based on the calculated metric.
#' @param field Character string specifying the data field to use for computing the metric.
#'        Acceptable values are `"abundance"`, `"rel_abundance"`, or `"norm_abundance"`.
#'        The field must not be empty in the `mgnet` object.
#' @param rank Optional character string specifying the taxonomic rank to consider when calculating metrics.
#'        If provided, the computation aggregates data at this taxonomic rank within the `lineage` data.
#'        Note: Aggregating data at different taxonomic levels may sum taxa classified together,
#'        potentially affecting the interpretability of normalized matrices.
#' @param order_fun Function used to compute the aggregate metric across the selected field. The default 
#'        is `sum`, but other functions like `mean` or `median` can be used to tailor the analysis.
#' @param decreasing Logical indicating if sorting should be in decreasing order. When `TRUE` (default),
#'        taxa with the highest metric values are considered top.
#'
#' @return For `mgnet` objects, a character vector of taxa IDs representing the top `n` taxa.
#'         For `mgnetList` objects, returns a list where each element is a character vector of top `n`
#'         samples for each respective `mgnet`.
#' @note Ensure that the specified field is not empty in the `mgnet` object before calling this function. 
#'       Failing to do so will result in an error. Aggregating data at different taxonomic levels
#'       is the result of summing up the same taxa together and not always it has sense for norm_abundance.
#' @export
#' @name top_n_taxa
#' @aliases top_n_taxa,mgnet-method top_n_taxa,mgnetList-method
setGeneric("top_n_taxa", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {standardGeneric("top_n_taxa")})

setMethod("top_n_taxa", "mgnet", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {
  
  if (!field %in% c("abundance", "rel_abundance", "norm_abundance")) {
    stop("Invalid field specified; choose from 'abundance', 'rel_abundance', 'norm_abundance'.")
  } else if(field == "abundance" & length(object@abundance) == 0) {
    stop("abundance slot cannot be empty if field is setted to 'abundance'")
  } else if(field == "rel_abundance" & length(object@rel_abundance) == 0) {
    stop("rel_bundance slot cannot be empty if field is setted to 'rel_abundance'")
  } else if(field == "norm_abundance" & length(object@norm_abundance) == 0) {
    stop("norm_abundance slot cannot be empty if field is setted to 'norm_abundance'")
  }
  
  if (!is.null(rank) && !rank %in% colnames(object@lineage)) {
    stop("Specified rank not found in lineage data.")
  } 
  
  # Extract the data matrix
  if(is.null(rank)){
    
    if(field == "abundance"){
      data_matrix <- abundance(object)
    } else if(field == "rel_abundance"){
      data_matrix <- rel_abundance(object)
    } else if(field == "norm_abundance"){
      data_matrix <- norm_abundance(object)
    }
    
  } else {
    
    if(field == "abundance"){
      data_matrix <- abundance(object, rank = rank)
    } else if(field == "rel_abundance"){
      data_matrix <- rel_abundance(object, rank = rank)
    } else if(field == "norm_abundance"){
      data_matrix <- norm_abundance(object, rank = rank)
    }
    
  }
  
  
  # Apply the order function to sum up or calculate the desired metric across columns (taxa)
  taxa_scores <- apply(data_matrix, 2, order_fun)
  
  # Order by the computed scores
  if (decreasing) {
    taxa_indices <- order(taxa_scores, decreasing = TRUE)
  } else {
    taxa_indices <- order(taxa_scores, decreasing = FALSE)
  }
  
  # Subset to top n
  top_n_indices <- taxa_indices[1:n]
  
  return(taxa_id(object)[top_n_indices])
  
})

setMethod("top_n_taxa", "mgnetList", function(object, n, field, rank = NULL, order_fun = sum, decreasing = TRUE) {
  
  result <- sapply(object, FUN = \(x){
    top_n_taxa(x, n = n, field = field, rank = rank, 
                  order_fun = order_fun, decreasing = decreasing)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(result)
})