#' Retrieve Abundance Data at Specified Taxonomic Rank
#'
#' Retrieves abundance data for an `mgnet` object, returning the abundance matrix
#' that represents the raw or processed counts of different taxa observed in each sample.
#' If a specific taxonomic rank is chosen, the function aggregates the abundance data 
#' by summing the counts of taxa classified under the same category at the specified rank,
#' providing a summarized view of abundance data at that rank. If no rank is specified,
#' or if the rank is missing, the function returns the original abundance matrix without 
#' aggregation. This allows for flexibility in analyzing data at various levels of 
#' taxonomic resolution.
#'
#' For an `mgnetList` object, it returns a list where each element corresponds to an 
#' `mgnet` object within the list, containing either the aggregated abundance data at 
#' the specified rank or the original abundance matrix if no rank is specified. This 
#' enables consistent data retrieval across multiple samples or experiments encapsulated 
#' within an `mgnetList`.
#'
#' This function differs from `taxa_id`, which refers to the unique identifiers 
#' associated with taxa at the finest taxonomic resolution (e.g., OTUs, species). 
#' While `taxa_id` often corresponds to the finest resolution, `taxa_name` at the same 
#' level might provide more descriptive names of the entities, and `abundance` at a 
#' specified rank leverages this taxonomic information to aggregate data accordingly.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param rank (Optional) A character string specifying the taxonomic rank of interest
#'        for data aggregation. If not provided or missing, the function returns the 
#'        original abundance data without aggregation.
#'
#' @return For an `mgnet` object, a numeric matrix of abundance data, aggregated at the 
#'         specified taxonomic rank if provided. For an `mgnetList` object, a list of numeric 
#'         matrices, each representing the aggregated or original abundance data from each 
#'         `mgnet` object within the list.
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang sym
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr select left_join group_by summarise
#' @name abundance
#' @aliases abundance,mgnet-method abundance,mgnetList-method
setGeneric("abundance", function(object, rank = "missing") standardGeneric("abundance"))

#' @rdname abundance
setMethod("abundance", signature(object = "mgnet", rank = "missing"), function(object) {
  object@abundance
})

#' @rdname abundance
setMethod("abundance", signature(object = "mgnet", rank = "character"), function(object, rank) {
  if(length(object@abundance) == 0 || !rank %in% colnames(object@lineage)) {
    stop("abundance and specified rank must be present in the object. Available ranks are: ", 
         paste(toString(colnames(object@lineage)), collapse=", "))
  }
  
  data_frame <- object@abundance %>% tibble::as_tibble(rownames="sample_id") %>%
    tidyr::pivot_longer(cols=-sample_id, names_to="taxa_id", values_to="abundance")
  
  lineage_info <- object@lineage %>%
    tibble::as_tibble(rownames="taxa_id") %>%
    dplyr::select(taxa_id, !!rlang::sym(rank))
  
  aggregated_data <- data_frame %>%
    dplyr::left_join(lineage_info, by = "taxa_id") %>%
    dplyr::group_by(sample_id, !!rlang::sym(rank)) %>%
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(rank), values_from = abundance)
  
  result_matrix <- aggregated_data %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  
  return(result_matrix)
})


#' @rdname abundance
setMethod("abundance", signature(object = "mgnetList", rank = "missing"), function(object) {
  sapply(object@mgnets, abundance, simplify=F, USE.NAMES=T)
})

#' @rdname abundance
setMethod("abundance", signature(object = "mgnetList", rank = "character"), function(object, rank) {
  sapply(object@mgnets, function(x) abundance(x, rank), simplify=F, USE.NAMES=T)
})
