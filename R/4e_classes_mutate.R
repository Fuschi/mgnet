# MUTATE INFO SAMPLE
#------------------------------------------------------------------------------#
#' Mutate Columns in info_sample of mgnet Objects
#'
#' This function allows users to modify existing columns or add new columns to the `info_sample` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's mutate semantics for flexible data transformation.
#'
#' @description
#' `mutate_info_sample` leverages dplyr's mutate functionality to enable column transformations or additions within the `info_sample`
#' data frame in `mgnet` objects.
#'
#' @param object An `mgnet` or `mgnetList` object to be modified.
#' @param ... Transformations to apply, passed to dplyr::mutate().
#'
#' @return The `mgnet` or `mgnetList` object with its `info_sample` slot updated.
#'
#' @seealso
#' \code{\link[dplyr]{mutate}} for details on transformation conditions.
#'
#' @export
#' @aliases mutate_info_sample,mgnet-method mutate_info_sample,mgnetList-method
#' @importFrom dplyr mutate
#' @importFrom rlang eval_tidy expr
#' @importFrom tibble column_to_rownames 
setGeneric("mutate_info_sample", function(object, ...) {standardGeneric("mutate_info_sample")})

setMethod("mutate_info_sample", "mgnet", function(object, ...) {
  
  # Check if info_sample is empty and initialize if necessary
  if(length(object@info_sample) == 0) {
    object@info_sample <- data.frame(sample_id = sample_id(object))
  }
  # Check for forbidden column name in the conditions
  new_names <- names(rlang::eval_tidy(rlang::expr(list(!!!conditions))))
  if("sample_id" %in% new_names) {
    stop("Modification or addition of 'sample_id' column is not allowed.")
  }
  
  conditions <- rlang::enquos(...)
  mutated_info_sample <- dplyr::mutate(info_sample(object, .fmt = "df"), !!!conditions)
  
  if("sample_id" %in% colnames(mutated_info_sample)){
    mutated_info_sample <- mutated_info_sample %>% column_to_rownames("sample_id")
  }
    
  object@info_sample <- mutated_info_sample 
  validObject(object)
  return(object)
})

setMethod("mutate_info_sample", "mgnetList", function(object, ...) {
  conditions <- rlang::enquos(...)
  object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
    
    if(length(mgnet_obj@info_sample) == 0) {
      mgnet_obj@info_sample <- data.frame(sample_id = sample_id(mgnet_obj))
    }
    
    mutated_info_sample <- dplyr::mutate(info_sample(mgnet_obj, .fmt = "df"), !!!conditions)
    
    if("sample_id" %in% colnames(mutated_info_sample)){
      mutated_info_sample <- mutated_info_sample %>% column_to_rownames("sample_id")
    }
    # Check for forbidden column name in the conditions
    new_names <- names(rlang::eval_tidy(rlang::expr(list(!!!conditions))))
    if("sample_id" %in% new_names) {
      stop("Modification or addition of 'sample_id' column is not allowed.")
    }
    
    mgnet_obj@info_sample <- mutated_info_sample
    validObject(mgnet_obj)
    return(mgnet_obj)
  })
  validObject(object)
  return(object)
})


# MUTATE INFO TAXA
#------------------------------------------------------------------------------#
#' Mutate Columns in info_taxa of mgnet Objects
#'
#' This function allows users to modify existing columns or add new columns to the `info_taxa` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's mutate semantics for flexible data transformation.
#'
#' @description
#' `mutate_info_taxa` leverages dplyr's mutate functionality to enable column transformations or additions within the `info_taxa`
#' data frame in `mgnet` objects. This can be particularly useful for adding new taxonomic information or adjusting existing ones.
#'
#' @param object An `mgnet` or `mgnetList` object to be modified.
#' @param ... Transformations to apply, passed to dplyr::mutate().
#' @return The `mgnet` or `mgnetList` object with its `info_taxa` slot updated.
#' @seealso
#' \code{\link[dplyr]{mutate}} for details on transformation conditions.
#' @export
#' @aliases mutate_info_taxa,mgnet-method mutate_info_taxa,mgnetList-method
#' @importFrom dplyr mutate
#' @importFrom rlang eval_tidy expr
#' @importFrom tibble column_to_rownames 
setGeneric("mutate_info_taxa", function(object, ...) {standardGeneric("mutate_info_taxa")})

setMethod("mutate_info_taxa", "mgnet", function(object, ...) {
  conditions <- rlang::enquos(...)
  
  # Check if info_taxa is empty and initialize if necessary
  if(length(object@info_taxa) == 0) {
    object@info_taxa <- data.frame(taxa_id = taxa_id(object))
  }
  # Check for forbidden column name in the conditions
  new_names <- names(rlang::eval_tidy(rlang::expr(list(!!!conditions))))
  if("taxa_id" %in% new_names) {
    stop("Modification or addition of 'taxa_id' column is not allowed.")
  }
  
  mutated_info_taxa <- dplyr::mutate(info_taxa(object, .fmt = "df"), !!!conditions)
  
  if("taxa_id" %in% colnames(mutated_info_taxa)){
    mutated_info_taxa <- mutated_info_taxa %>% column_to_rownames("taxa_id")
  }
  
  object@info_taxa <- mutated_info_taxa 
  validObject(object)
  return(object)
})

setMethod("mutate_info_taxa", "mgnetList", function(object, ...) {
  conditions <- rlang::enquos(...)
  object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
    
    # Check if info_taxa is empty and initialize if necessary
    if(length(mgnet_obj) == 0) {
      mgnet_obj@info_taxa <- data.frame(taxa_id = taxa_id(mgnet_obj))
    }
    # Check for forbidden column name in the conditions
    new_names <- names(rlang::eval_tidy(rlang::expr(list(!!!conditions))))
    if("taxa_id" %in% new_names) {
      stop("Modification or addition of 'taxa_id' column is not allowed.")
    }
    
    mutated_info_taxa <- dplyr::mutate(info_taxa(mgnet_obj, .fmt = "df"), !!!conditions)
    
    if("taxa_id" %in% colnames(mutated_info_taxa)){
      mutated_info_taxa <- mutated_info_taxa %>% column_to_rownames("taxa_id")
    }
    
    mgnet_obj@info_taxa <- mutated_info_taxa
    validObject(mgnet_obj)
    return(mgnet_obj)
  })
  validObject(object)
  return(object)
})
