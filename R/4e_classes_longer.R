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