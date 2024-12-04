# MGNET SUMMARY
#------------------------------------------------------------------------------#
#' Summarize `mgnet` and `mgnetList` Objects in Long Format
#'
#' Converts `mgnet` or `mgnetList` objects into a longer format, gathering all available 
#' information including abundances (absolute, relative, and normalized) and metadata into a 
#' single structured tibble. This transformation facilitates the use of tidyverse tools for
#' subsequent data manipulation and analysis.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        For `mgnet`, this method consolidates sample and taxa information along with 
#'        abundance data and associated metadata into a long-format tibble. For `mgnetList`,
#'        it applies the transformation to each `mgnet` object within the list and combines the results
#'        into a single tibble, incorporating an additional column named `mgnet` that identifies the source
#'        `mgnet` object.
#'
#' @details
#' The `gather_mgnet` method is tailored to prepare `mgnet` data for easy integration with tidyverse 
#' workflows by:
#' - Transforming data into a long format based on `sample_id` and `taxa_id`, making it more accessible for
#'   various tidyverse functions such as `mutate`, `filter`, and `summarize`.
#' - Aggregating both abundance metrics (e.g., `abun`, `rela`, `norm`) and metadata associated (e.g meta, taxa, comm_id) 
#'   with samples and taxa, which allows comprehensive data analysis within a single data frame structure.
#' - In the context of `mgnetList` objects, besides the standard data transformation applied to each `mgnet`,
#'   an additional column using the reserved keyword `mgnet` is included to keep track of the original `mgnet` 
#'   object each row of data pertains to, facilitating analysis across multiple datasets.
#'
#' @return A tibble in a long format that consolidates all the relevant data fields under sample and taxon identifiers.
#'         For `mgnetList` objects, the tibble includes an extra column `mgnet` that indicates the source `mgnet` object.
#'
#' @export
#' @aliases gather_mgnet,mgnet-method gather_mgnet,mgnetList-method
#' @importFrom tidyr pivot_longer expand_grid
#' @importFrom dplyr left_join mutate
#' @importFrom purrr map imap list_rbind
setGeneric("gather_mgnet", function(object) standardGeneric("gather_mgnet"))


setMethod("gather_mgnet", "mgnet", function(object){
  
  # empty
  if(miss_sample(object) || miss_taxa(object)) stop("Error: No sample or taxa available.")
  

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

setMethod("gather_mgnet", "mgnetList", function(object) {
  
  # empty
  if( miss_sample(object, "any") || miss_taxa(object, "any") ) stop("Error: No sample or taxa available in at least one of the mgnet objects.")
  
  purrr::map(object, gather_mgnet) %>%
    purrr::imap(\(x,y){
      x <- dplyr::mutate(x, mgnet = y, .before = 1)
    }) %>%
    purrr::list_rbind() %>%
    return()
  
})


#' Summarize Sample Metadata in `mgnet` and `mgnetList` Objects
#'
#' Converts `mgnet` or `mgnetList` objects into a longer format, focusing exclusively on 
#' sample metadata. This transformation facilitates the use of tidyverse tools for
#' subsequent data manipulation and analysis, especially when only sample-specific details are required.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        For `mgnet`, this method consolidates all sample metadata into a long-format tibble. 
#'        For `mgnetList`, it applies the transformation to each `mgnet` object within the list and combines the results
#'        into a single tibble, incorporating an additional column named `mgnet` that identifies the source
#'        `mgnet` object.
#'
#' @details
#' The `gather_meta` method is tailored to prepare `mgnet` data for easy integration with tidyverse 
#' workflows by:
#' - Transforming sample metadata into a long format based on `sample_id`, making it more accessible for
#'   various tidyverse functions such as `mutate`, `filter`, and `summarize`.
#' - For `mgnetList` objects: The transformation is applied to each individual `mgnet` in the list. 
#'   The resulting longer format data for each `mgnet` is then combined into a single tibble, with 
#'   an additional column using the reserved keyword `mgnet` to track the original `mgnet` object each row of data pertains to.
#'
#' @return A tibble in a long format that consolidates all the relevant sample metadata under sample identifiers.
#'         For `mgnetList` objects, the tibble includes an extra column `mgnet` that indicates the source `mgnet` object.
#'
#' @export
#' @aliases gather_meta,mgnet-method gather_meta,mgnetList-method
#' @importFrom dplyr left_join mutate
#' @importFrom purrr map imap list_rbind
setGeneric("gather_meta", function(object) standardGeneric("gather_meta"))

setMethod("gather_meta", "mgnet", function(object){
  
  if(miss_sample(object)) stop("Error: No sample available.")
  
  if (has_slot(object, "meta")) {
    return(meta(object, .fmt = "tbl"))
  } else {
    return(tibble::tibble(sample_id = sample_id(object)))
  }
  
})

setMethod("gather_meta", "mgnetList", function(object) {
  
  if(miss_sample(object, "any")) stop("Error: No sample available in at least one of the mgnet objects.")
  
  object %>%
    purrr::map(gather_meta) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind() %>%
    return()
  
})


#' Summarize Taxa Metadata in `mgnet` and `mgnetList` Objects
#'
#' Converts `mgnet` or `mgnetList` objects into a longer format, focusing exclusively on 
#' taxa metadata. This transformation facilitates the use of tidyverse tools for
#' subsequent data manipulation and analysis, especially when only sample-specific details are required.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        For `mgnet`, this method consolidates all taxa metadata into a long-format tibble. 
#'        For `mgnetList`, it applies the transformation to each `mgnet` object within the list and combines the results
#'        into a single tibble, incorporating an additional column named `mgnet` that identifies the source
#'        `mgnet` object.
#'
#' @details
#' The `gather_taxa` method is tailored to prepare `mgnet` data for easy integration with tidyverse 
#' workflows by:
#' - Transforming taxa metadata into a long format based on `taxa_id`, making it more accessible for
#'   various tidyverse functions such as `mutate`, `filter`, and `summarize`.
#' - For `mgnetList` objects: The transformation is applied to each individual `mgnet` in the list. 
#'   The resulting longer format data for each `mgnet` is then combined into a single tibble, with 
#'   an additional column using the reserved keyword `mgnet` to track the original `mgnet` object each row of data pertains to.
#'
#' @return A tibble in a long format that consolidates all the relevant taxa metadata under sample identifiers.
#'         For `mgnetList` objects, the tibble includes an extra column `mgnet` that indicates the source `mgnet` object.
#'
#' @export
#' @aliases gather_taxa,mgnet-method gather_taxa,mgnetList-method
#' @importFrom dplyr left_join mutate
#' @importFrom purrr map imap list_rbind
setGeneric("gather_taxa", function(object) standardGeneric("gather_taxa"))

setMethod("gather_taxa", "mgnet", function(object){
  
  if(miss_taxa(object)) stop("Error: No taxa available.")
  
  if (has_metataxa(object)) {
    return(taxa(object, .fmt = "tbl"))
  } else {
    return(tibble::tibble(taxa_id = taxa_id(object)))
  }
  
})

setMethod("gather_taxa", "mgnetList", function(object) {
  
  if(miss_taxa(object, "any")) stop("Error: No taxa available in at least one of the mgnet objects.")
  
  object %>%
    purrr::map(gather_taxa) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind() %>%
    return()
  
})


#' Retrieve Edge List with Metadata from mgnet Object(s)
#'
#' Extracts an edge list from the network slot of an `mgnet` or `mgnetList` object, including edge attributes 
#' and associated taxa metadata for both vertices in each link.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'               For `mgnet`, the method will extract the network's edge list and merge taxa metadata 
#'               for both source and target nodes.
#'               For `mgnetList`, the method will apply the transformation to each `mgnet` in the list 
#'               and bind the results, including an additional column `mgnet` to track the source object.
#' @param .suffix A character vector of length 2 providing suffixes to append to the taxa 
#'        metadata columns to distinguish the nodes connected from an edge. 
#'        Default is c("_from", "_to"). The values must be distinct to prevent column name overlap.
#'
#' @details
#' This method leverages the network structure within `mgnet` objects to generate a detailed edge tibble that includes:
#' - `Node Identifiers`: Each identifier for the nodes connected by an edge begins with the string "taxa" and is 
#'   appended with `.suffix` to indicate the source and target nodes, respectively (e.g., `taxa_id_1`, `taxa_id_2`).
#' - `Link Attributes`: Includes all available attributes associated with the links, such as weight.
#' - `Metadata`: Metadata from both source and target taxa are included, with column names appended 
#'   with `_1` and `_2` as suffixes. These suffixes help distinguish between the source and target taxa metadata.
#'
#' For `mgnetList` objects, the transformation is applied individually to each `mgnet` object, and the results 
#' are combined into a single dataframe with an additional column named `mgnet` that identifies the source object.
#'
#' @return A dataframe containing the edge list with additional metadata columns for both vertices involved in each link.
#'         For `mgnetList` objects, the output is a single dataframe where results from all `mgnet` objects are bound together.
#'
#' @importFrom igraph E is_weighted as_data_frame
#' @importFrom dplyr left_join
#' @importFrom purrr map imap list_rbind
#' @export
#' @name gather_link
#' @aliases gather_link,mgnet-method gather_link,mgnetList-method
setGeneric("gather_link", function(object, .suffix = c("_from", "_to")) standardGeneric("gather_link"))

setMethod("gather_link", "mgnet", function(object, .suffix = c("_from", "_to")) {
  
  link(object)
  
})

setMethod("gather_link", "mgnetList", function(object, .suffix = c("_from", "_to")) {
  
  purrr::map(object, link, .suffix = .suffix) %>% 
    purrr::imap(~ dplyr::mutate(.x, mgnet = .y, .before = 1)) %>%  
    purrr::list_rbind() 
  
})

