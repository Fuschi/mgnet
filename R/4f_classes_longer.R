#------------------------------------------------------------------------------#
#' Transform mgnet Object(s) to Longer Format
#'
#' Converts an `mgnet` or `mgnetList` object into a "longer" format, collecting all possible 
#' information (abundance, relative abundance, normalized abundance, and metadata) into a single 
#' tibble. This method uses a tidyverse approach to gather data into a structure suitable for 
#' further analysis.
#'
#' @param object An `mgnet` or `mgnetList` object. For `mgnet`, the method will transform the data 
#'        to a longer format by combining information on samples, taxa, abundances, 
#'        and metadata. For `mgnetList`, the method will apply the transformation to each `mgnet` 
#'        object and bind the results into a single tibble, with an additional column named `mgnet` 
#'        indicating the source `mgnet` object.
#'
#' @details 
#' This method utilizes specific reserved keywords from the `mgnet` class, including:
#' - `sample_id`: The identifier for each sample.
#' - `taxa_id`: The identifier for each taxon.
#' - `comm_id`: The identifier for membership community for each taxon.
#' - `mgnet`: Added when dealing with `mgnetList`, to track the source `mgnet` object.
#' - `abun`: Refers to abundance values in the data.
#' - `rela`: Refers to relative abundance values in the data.
#' - `norm`: Refers to normalized abundance values in the data.
#' 
#' - For `mgnet` objects: The method gathers data on samples (`sample_id`), taxa (`taxa_id`), 
#'   and abundance values (e.g., `abun`, `rela`, `norm`) into a longer format, combining sample 
#'   metadata and taxa metadata as well.
#' - For `mgnetList` objects: The transformation is applied to each individual `mgnet` in the list. 
#'   The resulting longer format data for each `mgnet` is then combined into a single tibble, with 
#'   an additional column `mgnet`, which contains the name of each `mgnet` object. This allows 
#'   users to easily track which data comes from which `mgnet`.
#' - If the network community information are available a column `comm_id` is 
#'   added with the membership for each taxon.
#'
#' The tidyverse functions used include `pivot_longer` for reshaping data, and `left_join` for 
#' merging sample and taxa information. The `mgnet_longer` method ensures that all relevant data 
#' (e.g., abundance, relative abundance, normalized abundance, metadata) is organized in a 
#' consistent and accessible manner.
#'
#' @return A tibble in a longer format that includes all the combined information. For `mgnetList` 
#'         objects, the output is a single tibble where results from all `mgnet` objects are bound 
#'         together, and an extra `mgnet` column identifies the source `mgnet` object.
#'         
#' @importFrom tidyr pivot_longer expand_grid
#' @importFrom dplyr left_join mutate join_by
#' @importFrom purrr map imap list_rbind
#' @export
#' @aliases mgnet_longer,mgnet-method mgnet_longer,mgnetList-method
setGeneric("mgnet_longer", function(object) standardGeneric("mgnet_longer"))


setMethod("mgnet_longer", "mgnet", function(object){
  
  # empty
  if( nsample(object)==0 | ntaxa(object)==0 ) stop("object must have samples and taxa")
  

  long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                     taxa_id = taxa_id(object))
  if(length(object@abun)!=0) {
   long_mgnet <- long_mgnet %>% 
     dplyr::left_join(
       abun(object, .fmt="tbl") %>%
         tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abun"),
         dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(length(object@rela)!=0) {
   long_mgnet <- long_mgnet %>% 
     dplyr::left_join(
       rela(object, .fmt="tbl") %>%
         tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "rela"),
       dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(length(object@norm)!=0) {
   long_mgnet <- long_mgnet %>% 
     dplyr::left_join(
       norm(object, .fmt="tbl") %>%
         tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "norm"),
       dplyr::join_by("sample_id", "taxa_id"))
  }
  
  if(length(meta(object))!=0){
   long_mgnet <- long_mgnet %>% 
     dplyr::left_join(
       meta(object, .fmt="tbl"), by = "sample_id")
  }
  
  if(length(taxa(object))!=0){
   long_mgnet <- long_mgnet %>% 
     dplyr::left_join(
       taxa(object, .fmt="tbl"), by = "taxa_id")
  }
  
  return(long_mgnet)
  
})

setMethod("mgnet_longer", "mgnetList", function(object) {
  
  # empty
  if( any(nsample(object))==0 | any(ntaxa(object))==0 ) stop("all mgnet in object must have samples and taxa")
  
  purrr::map(object, mgnet_longer) %>%
    purrr::imap(\(x,y){
      x <- dplyr::mutate(x, mgnet = y, .before = 1)
    }) %>%
    purrr::list_rbind() %>%
    return()
  
})


#' Retrieve Edge List with Metadata from mgnet Object(s)
#'
#' Extracts an edge list from the network slot of an `mgnet` or `mgnetList` object, including edge weights 
#' and associated taxa metadata for both vertices in each link.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'               For `mgnet`, the method will extract the network's edge list and merge taxa metadata 
#'               for both source and target nodes.
#'               For `mgnetList`, the method will apply the transformation to each `mgnet` in the list 
#'               and bind the results, including an additional column `mgnet` to track the source object.
#'
#' @details
#' This method utilizes the network structure within `mgnet` to generate a comprehensive edge list 
#' that includes:
#' - `source`: The identifier for the source node.
#' - `target`: The identifier for the target node.
#' - `weight`: The weight of the edge, if available.
#' - Metadata from both source and target taxa, added as additional columns.
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
#' @name mgnet_edge_list
#' @aliases mgnet_edge_list,mgnet-method mgnet_edge_list,mgnetList-method
setGeneric("mgnet_edge_list", function(object) standardGeneric("mgnet_edge_list"))

setMethod("mgnet_edge_list", "mgnet", function(object) {
  
  # Ensure the network is available
  if (is.null(object@netw)) {
    stop("The mgnet object does not contain a network.")
  }
  
  # Extract edge list with weights
  net <- netw(object)
  edges_df <- igraph::as_data_frame(net, what = "edges")
  
  # If the network has weights, include them
  if ("weight" %in% igraph::is_weighted(net)) {
    edges_df$weight <- igraph::E(net)$weight
  }
  
  # Check if the network is directed
  directed <- igraph::is.directed(net)
  
  # Merge with taxa metadata
  if(length(taxa(object)) != 0 & directed){
    
    taxa_info <- taxa(object, .fmt = "tbl")
    colnames(taxa_info) <- paste0(colnames(taxa_info), "_from")
    edges_df <- edges_df %>%
      left_join(taxa_info, by = c("from" = "taxa_id_from")) 
    
    taxa_info <- taxa(object, .fmt = "tbl")
    colnames(taxa_info) <- paste0(colnames(taxa_info), "_to")
    edges_df <- edges_df %>%
      left_join(taxa_info, by = c("from" = "taxa_id_to")) 
    
  } else if(length(taxa(object)) != 0 & !directed){
    
    colnames(edges_df)[1:2] <- c("node1", "node2")
    
    taxa_info <- taxa(object, .fmt = "tbl")
    colnames(taxa_info) <- paste0(colnames(taxa_info), "_1")
    edges_df <- edges_df %>%
      left_join(taxa_info, by = c("node1" = "taxa_id_1")) 
    
    taxa_info <- taxa(object, .fmt = "tbl")
    colnames(taxa_info) <- paste0(colnames(taxa_info), "_2")
    edges_df <- edges_df %>%
      left_join(taxa_info, by = c("node2" = "taxa_id_2")) 
  }
  
  return(edges_df)
})

setMethod("mgnet_edge_list", "mgnetList", function(object) {
  
  # Apply edge list extraction to each mgnet in the list
  purrr::map(object, mgnet_edge_list) %>%
    purrr::imap(\(x, y) {
      x <- dplyr::mutate(x, mgnet = y, .before = 1)
    }) %>%
    purrr::list_rbind() 
  
})