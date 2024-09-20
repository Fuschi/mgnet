# TOP_SAMPLE
#------------------------------------------------------------------------------#
#' Return Ordered the Sample IDs Based on Aggregate Metrics on Abundances Matrices
#'
#' Return the ordered samples from an `mgnet` object based on an aggregate metric calculated from
#' specified data fields such as abundance, relative abundance, or normalized abundance, with
#' the chosen metric, which can be customized.
#'
#' @param object An `mgnet` object.
#' @param field Character string specifying the data field to use for computing the metric.
#'        Acceptable values are `"abun"`, `"rela"`, or `"norm"`.
#'        The field must not be empty in the `mgnet` object.
#' @param .var A character string specifying the column name from the `taxa` slot to be used
#'        for grouping and aggregation. If this parameter is omitted, the function will return
#'        the original abundance data without any aggregation.
#' @param .aggregation_fun A function to specify how the abundance data should be aggregated.
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#' @param .order_fun Function used to compute the aggregate metric across the selected field. The default
#'        is `sum`, but other functions like `mean` or `median` can be used to tailor the analysis.
#' @param decreasing Logical indicating if sorting should be in decreasing order. When `TRUE` (default),
#'        samples with the highest metric values are considered top.
#'
#' @return For `mgnet` objects, a character vector of sample IDs representing the ordered samples.
#'         For `mgnetList` objects, returns a list where each element is a character vector of
#'         the ordered samples for each respective `mgnet`.
#' @note Ensure that the specified field is not empty in the `mgnet` object before calling this function.
#'       Failing to do so will result in an error.
#' @export
#' @name top_samples
#' @aliases top_samples,mgnet-method top_samples,mgnetList-method
setGeneric("top_samples", function(object, field, .var, 
                                   .aggregation_fun = sum, .order_fun = sum, decreasing = TRUE) {standardGeneric("top_samples")})

setMethod("top_samples", "mgnet", function(object, field, .var,
                                           .aggregation_fun = sum, .order_fun = sum, decreasing = TRUE) {

  if (!field %in% c("abun", "rela", "norm")) {
    stop("Invalid field specified; choose from 'abun', 'rela', 'norm'.")
  } else if(field == "abun" & length(object@abun) == 0) {
    stop("abun slot cannot be empty if field is setted to 'abun'")
  } else if(field == "rela" & length(object@rela) == 0) {
    stop("rel_bundance slot cannot be empty if field is setted to 'rela'")
  } else if(field == "norm" & length(object@norm) == 0) {
    stop("norm slot cannot be empty if field is setted to 'norm'")
  }

  if(!missing(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name in taxa.")
  }

  # Extract the data matrix
  if(missing(.var)){

    if(field == "abun"){
      data_matrix <- abun(object)
    } else if(field == "rela"){
      data_matrix <- rela(object)
    } else if(field == "norm"){
      data_matrix <- norm(object)
    }

  } else {

    if(field == "abun"){
      data_matrix <- abun(object, .var = .var, .fun = .aggregation_fun)
    } else if(field == "rela"){
      data_matrix <- rela(object, .var = .var, .fun = .aggregation_fun)
    } else if(field == "norm"){
      data_matrix <- norm(object, .var = .var, .fun = .aggregation_fun)
    }

  }

  # Apply the order function to sum up or calculate the desired metric across columns (taxa)
  sample_scores <- apply(data_matrix, 1, .order_fun)

  # Order by the computed scores
  if (decreasing) {
    sample_indices <- order(sample_scores, decreasing = TRUE)
  } else {
    sample_indices <- order(sample_scores, decreasing = FALSE)
  }


  return(rownames(data_matrix)[sample_indices])
})

setMethod("top_samples", "mgnetList", function(object, field, .var, 
                                               .aggregation_fun = sum, .order_fun = sum, decreasing = TRUE) {

  result <- sapply(object, FUN = \(x){
    top_samples(x, field = field, .var = .var, .aggregation_fun = .aggregation_fun,
                  .order_fun = .order_fun, decreasing = decreasing)
  }, simplify = FALSE, USE.NAMES = TRUE)

  return(result)
})


# TOP_TAXA
#------------------------------------------------------------------------------#
#' Return Ordered the Taxa IDs Based on Aggregate Metrics on Abundances Matrices
#'
#' Return the ordered taxa from an `mgnet` object based on an aggregate metric calculated from
#' specified data fields such as abundance, relative abundance, or normalized abundance, with
#' the chosen metric, which can be customized.
#'
#' @param object An `mgnet` object.
#' @param field Character string specifying the data field to use for computing the metric.
#'        Acceptable values are `"abun"`, `"rela"`, or `"norm"`.
#'        The field must not be empty in the `mgnet` object.
#' @param .var A character string specifying the column name from the `taxa` slot to be used
#'        for grouping and aggregation. If this parameter is omitted, the function will return
#'        the original abundance data without any aggregation.
#' @param .aggregation_fun A function to specify how the abundance data should be aggregated.
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#' @param .order_fun Function used to compute the aggregate metric across the selected field. The default
#'        is `sum`, but other functions like `mean` or `median` can be used to tailor the analysis.
#' @param decreasing Logical indicating if sorting should be in decreasing order. When `TRUE` (default),
#'        taxa with the highest metric values are considered top.
#'
#' @return For `mgnet` objects, a character vector of sample IDs representing the ordered taxa.
#'         For `mgnetList` objects, returns a list where each element is a character vector of
#'         the ordered taxa for each respective `mgnet`.
#' @note Ensure that the specified field is not empty in the `mgnet` object before calling this function.
#'       Failing to do so will result in an error.
#' @export
#' @name top_taxa
#' @aliases top_taxa,mgnet-method top_taxa,mgnetList-method
setGeneric("top_taxa", function(object, field, .var, .aggregation_fun = sum,.order_fun = sum, decreasing = TRUE) {standardGeneric("top_taxa")})

setMethod("top_taxa", "mgnet", function(object, field, .var, .aggregation_fun = sum, .order_fun = sum, decreasing = TRUE) {
  
  if (!field %in% c("abun", "rela", "norm")) {
    stop("Invalid field specified; choose from 'abun', 'rela', 'norm'.")
  } else if(field == "abun" & length(object@abun) == 0) {
    stop("abun slot cannot be empty if field is setted to 'abun'")
  } else if(field == "rela" & length(object@rela) == 0) {
    stop("rel_bundance slot cannot be empty if field is setted to 'rela'")
  } else if(field == "norm" & length(object@norm) == 0) {
    stop("norm slot cannot be empty if field is setted to 'norm'")
  }
  
  if(!missing(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  # Extract the data matrix
  if(missing(.var)){
    
    if(field == "abun"){
      data_matrix <- abun(object)
    } else if(field == "rela"){
      data_matrix <- rela(object)
    } else if(field == "norm"){
      data_matrix <- norm(object)
    }
    
  } else {
    
    if(field == "abun"){
      data_matrix <- abun(object, .var = .var, .fun = .aggregation_fun)
    } else if(field == "rela"){
      data_matrix <- rela(object, .var = .var, .fun = .aggregation_fun)
    } else if(field == "norm"){
      data_matrix <- norm(object, .var = .var, .fun = .aggregation_fun)
    }
    
  }
  
  # Apply the order function to sum up or calculate the desired metric across columns (taxa)
  taxa_scores <- apply(data_matrix, 2, .order_fun)
  
  # Order by the computed scores
  if (decreasing) {
    taxa_indices <- order(taxa_scores, decreasing = TRUE)
  } else {
    taxa_indices <- order(taxa_indices, decreasing = FALSE)
  }
  
  
  return(colnames(data_matrix)[taxa_indices])
})

setMethod("top_taxa", "mgnetList", function(object, field, .var, .aggregation_fun, .order_fun = sum, decreasing = TRUE) {
  
  result <- sapply(object, FUN = \(x){
    top_taxa(x, field = field, .var = .var, .aggregation_fun = .aggregation_fun,
                .order_fun = .order_fun, decreasing = decreasing)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(result)
})