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


#' @title Retrieve Edge List with Metadata from mgnet Objects
#'
#' @description
#' Extracts edge information from the network slot of an `mgnet` or `mgnetList` object,
#' optionally merging metadata from both vertices (nodes) into each link. 
#' Internally uses \code{igraph::as_data_frame()} to get the edge list, then
#' joins with taxa metadata if \code{.suffix} is provided.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'   For `mgnet`, the method extracts the network's edges and merges taxa metadata
#'   (if requested).
#'   For `mgnetList`, the method applies the transformation to each `mgnet` in
#'   the list and combines the results, adding a `mgnet` column for the source.
#' @param .suffix A character vector of length 2 providing suffixes to append
#'   to the taxa metadata columns for the two nodes in each link. The two suffixes
#'   must be distinct (e.g. \code{c("_from", "_to")}). If missing or \code{NULL},
#'   no metadata join is performed.
#'
#' @details
#' - **Edge Attributes**: All available attributes on the edges (e.g., weight) are
#'   preserved.  
#' - **Node Metadata**: If \code{.suffix} is given, columns from the \code{from} and 
#'   \code{to} nodes' metadata are appended with the respective suffixes, so they
#'   do not overwrite each other or edge attributes.  
#' - **Filtering**: If the \code{mgnet} object has selected links (via \code{select_link}),
#'   only those links are kept.
#'
#' @return A dataframe (tibble) containing the edges plus node-level metadata
#'   (if requested). For \code{mgnetList}, the result is a single tibble with
#'   an extra \code{mgnet} column identifying each sub-object.
#'
#' @importFrom igraph E as_data_frame
#' @importFrom dplyr left_join mutate
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble as_tibble
#' @export
#' @name gather_link
#' @aliases gather_link,mgnet-method gather_link,mgnetList-method
setGeneric("gather_link", function(object, .suffix) standardGeneric("gather_link"))


#----------------#
# mgnet Method  #
#----------------#
setMethod("gather_link", "mgnet", function(object, .suffix) {
  
  ## 1) Check that the network slot is available
  if (miss_slot(object, "netw")) {
    stop("Error: No network available in this mgnet object.")
  }
  
  ## 2) Subset edges if there's a selection
  selected_links <- get_selected_links(object)
  if (!is.null(selected_links)) {
    # create subgraph from selected edges only
    netw(object) <- igraph::subgraph_from_edges(
      graph = netw(object),
      eids  = selected_links,
      delete.vertices = FALSE
    )
  }
  
  ## 3) Extract the edges using igraph::as_data_frame()
  net      <- netw(object)
  edges_df <- tibble::as_tibble(igraph::as_data_frame(net, what = "edges"))
  
  ## 4) If .suffix is provided, do checks and merge node metadata
  if (!missing(.suffix) && !is.null(.suffix)) {
    
    # Check .suffix length, type, distinctness
    if (length(.suffix) != 2 || !is.character(.suffix) || .suffix[1] == .suffix[2]) {
      stop(".suffix must be a character vector of length 2 with distinct values.")
    }
    
    # Check for overlap between suffix-based node columns and existing edge columns
    # We haven't created them yet, but we'll see if they'd conflict with edge columns:
    # e.g., "family_from", "family_to" might already exist in edges_df.
    # 1) Build the potential new column names for from/to node metadata
    from_cols <- paste0(taxa_vars(object), .suffix[1])
    to_cols   <- paste0(taxa_vars(object), .suffix[2])
    new_taxa_info_cols <- c(from_cols, to_cols)
    
    # 2) If any of these new columns already appear in edges_df, that's a conflict
    conflict_cols <- intersect(new_taxa_info_cols, names(edges_df))
    if (length(conflict_cols) > 0) {
      stop(
        "Column name conflict: the following columns already exist in edges_df: ",
        paste(conflict_cols, collapse = ", "),
        ". Please pick different suffixes or rename edge columns."
      )
    }
    
    # 3) Also check that no existing edge attributes end with those suffixes
    #    (if you want to enforce that rule). Only relevant if you specifically
    #    want to disallow edge attributes named "weight_from", for instance.
    edges_attr_names <- names(igraph::edge_attr(net))
    if (!is.null(edges_attr_names) && length(edges_attr_names) > 0) {
      # Build a small pattern like:  (?:_from|_to)$
      # This pattern checks for any string that ends with one of the suffixes
      pattern <- paste0("(", paste0(.suffix, collapse="|"), ")$")
      invalid_names <- edges_attr_names[stringr::str_detect(edges_attr_names, pattern)]
      if (length(invalid_names) > 0) {
        stop(
          "Edge attributes end with disallowed suffixes (", paste(.suffix, collapse=", "), 
          "): ", paste(invalid_names, collapse=", "),
          ". Suffixes must not match existing edge attribute endings."
        )
      }
    }
    
    # 4) Merge with node metadata for source nodes
    tbl_from <- gather_taxa(object)
    colnames(tbl_from)[-1] <- paste0(colnames(tbl_from)[-1], .suffix[1])
    edges_df <- dplyr::left_join(edges_df, tbl_from, 
                                 by = dplyr::join_by(from == taxa_id))
    
    # 5) Merge with node metadata for target nodes
    tbl_to <- gather_taxa(object)
    colnames(tbl_to)[-1] <- paste0(colnames(tbl_to), .suffix[2])[-1]
    edges_df <- dplyr::left_join(edges_df, tbl_to, 
                                 by = dplyr::join_by(to == taxa_id))
  }
  
  ## 5) Return the final edges df
  edges_df
})


#------------------#
# mgnetList Method #
#------------------#
setMethod("gather_link", "mgnetList", function(object, .suffix) {
  
  # For each mgnet in the list, gather_link() individually
  # .suffix might be missing or provided
  result <- purrr::imap(object, ~ {
    if (missing(.suffix) || is.null(.suffix)) {
      gather_link(.x)  # no suffix
    } else {
      gather_link(.x, .suffix = .suffix)
    }
  })
  
  # Add a mgnet column to each tibble, then bind them
  # .x is the tibble, .y is the mgnet name
  combined <- purrr::imap(result, ~ dplyr::mutate(.x, mgnet = .y, .before = 1)) %>%
    purrr::list_rbind()
  
  combined
})


