#' Summarize `mgnet` and `mgnets` Objects in Long Format
#'
#' Converts `mgnet` or `mgnets` objects into a longer format including abundances 
#' (absolute, relative, and normalized) and metadata into a single structured tibble. 
#' This transformation facilitates the use of tidyverse tools for
#' subsequent data manipulation and analysis.
#'
#' @param object An `mgnet` or `mgnets` object.
#'        For `mgnet`, this method consolidates sample and taxa information along with 
#'        abundance data and associated metadata into a long-format tibble. For `mgnets`,
#'        it applies the transformation to each `mgnet` object within the list and combines the results
#'        into a single tibble, incorporating an additional column named `mgnet` that identifies the source
#'        `mgnet` object.
#'
#' @details
#' The `mgnet_longer` method is tailored to prepare `mgnet` data for easy integration with tidyverse 
#' workflows by:
#' - Transforming data into a long format based on `sample_id` and `taxa_id`, making it more accessible for
#'   various tidyverse functions such as `mutate`, `filter`, and `summarize`.
#' - Aggregating both abundance metrics (e.g., `abun`, `rela`, `norm`) and metadata associated (e.g meta, taxa, comm_id) 
#'   with samples and taxa, which allows comprehensive data analysis within a single data frame structure.
#' - In the context of `mgnets` objects, besides the standard data transformation applied to each `mgnet`,
#'   an additional column using the reserved keyword `mgnet` is included to keep track of the original `mgnet` 
#'   object each row of data pertains to, facilitating analysis across multiple datasets.
#'
#' @return A tibble in a long format that consolidates all the relevant data fields under sample and taxon identifiers.
#'         For `mgnets` objects, the tibble includes an extra column `mgnet` that indicates the source `mgnet` object.
#'
#' @export
#' @aliases mgnet_longer,mgnet-method mgnet_longer,mgnets-method
setGeneric("mgnet_longer", function(object) standardGeneric("mgnet_longer"))


setMethod("mgnet_longer", "mgnet", function(object){
  
  # empty
  if(miss_sample(object) || miss_taxa(object)) cli::cli_abort("Error: No sample or taxa available.")
  
  
  long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                   taxa_id = taxa_id(object))
  if(has_slot(object, "abun")) {
    long_mgnet <- long_mgnet %>% 
      dplyr::left_join(
        abun(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abun"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "rela")) {
    long_mgnet <- long_mgnet %>% 
      dplyr::left_join(
        rela(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "rela"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "norm")) {
    long_mgnet <- long_mgnet %>% 
      dplyr::left_join(
        norm(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "norm"),
        dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(has_slot(object, "meta")){
    long_mgnet <- long_mgnet %>% 
      dplyr::left_join(
        meta(object, .fmt="tbl"), by = "sample_id")
  }
  
  if(has_slot(object, "taxa")){
    long_mgnet <- long_mgnet %>% 
      dplyr::left_join(
        taxa(object, .fmt="tbl"), by = "taxa_id")
  }
  
  return(long_mgnet)
  
})

setMethod("mgnet_longer", "mgnets", function(object) {
  
  # empty
  if( miss_sample(object, "any") || miss_taxa(object, "any") ) cli::cli_abort("Error: No sample or taxa available in at least one of the mgnet objects.")
  
  purrr::map(object, mgnet_longer) %>%
    purrr::imap(\(x,y){
      x <- dplyr::mutate(x, mgnet = y, .before = 1)
    }) %>%
    purrr::list_rbind() %>%
    return()
  
})