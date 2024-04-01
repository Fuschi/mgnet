#' Convert mgnetList Object to Longer Format
#'
#' This function transforms data from each `mgnet` object within an `mgnetList` into a longer format,
#' combining abundance data, sample metadata, taxonomic information, and community memberships into a
#' single, comprehensive tibble. Each row represents a unique combination of sample and taxa, enriched
#' with the specified data types. For `mgnetList` objects, an additional column is added to identify
#' the source `mgnet` object for each row, facilitating data origin tracking. The `source_mgnet_to`
#' parameter is mandatory for `mgnetList` objects to specify this column name.
#'
#' @param object An `mgnetList` object containing multiple `mgnet` objects to be combined and
#'        transformed into a longer format.
#' @param source_mgnet_to A string specifying the column name to be added for identifying the source 
#'        `mgnet` object within an `mgnetList`. This parameter is required for `mgnetList` objects and
#'        must be provided; it is not used for single `mgnet` objects.
#' @param abundance A logical value; if `TRUE`, includes abundance data for each taxa in each sample.
#' @param rel_abundance A logical value; if `TRUE`, includes relative abundance data for each taxa
#'        in each sample, facilitating comparison across samples.
#' @param norm_abundance A logical value; if `TRUE`, includes normalized abundance data for each taxa,
#'        useful for downstream analyses that require standardized input.
#' @param lineage A logical value; if `TRUE`, includes lineage information for each taxa, enriching
#'        the dataset with taxonomic context.
#' @param info_sample A logical value; if `TRUE`, includes sample metadata, such as environmental
#'        conditions or collection details, enhancing interpretability.
#' @param info_taxa A logical value; if `TRUE`, includes detailed taxonomic information for each taxa,
#'        supporting fine-grained biological insights.
#' @param community A logical value; if `TRUE`, includes community membership data for each taxa, 
#'        enabling analysis of community structure and diversity.
#'
#' @return For `mgnetList` objects, returns a `tibble` in a long format containing combined data from 
#'         each `mgnet` object in the list, enriched with specified metadata and abundance information,
#'         and an additional column indicating the source `mgnet` object as specified by `source_mgnet_to`.
#'         For `mgnet` objects, returns a similar `tibble` without the source identification column.
#'
#' @export
#' @aliases mgnet_longer,mgnetList-method
#'
#' @importFrom tibble tibble enframe
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select right_join bind_rows mutate 
#' @importFrom rlang enquos :=
#'
#' @note The function checks for the presence of requested data types in each `mgnet` object and
#'       automatically excludes any data types not present or specifically set to `FALSE`.
#'       It's essential to ensure that each `mgnet` object within the `mgnetList` is properly
#'       formatted and contains the necessary data for the requested outputs.
setGeneric("mgnet_longer", function(object, source_mgnet_to,
                                    abundance = NULL, rel_abundance = NULL, norm_abundance = NULL,
                                    lineage = NULL, info_sample = NULL, info_taxa = NULL, community = NULL) standardGeneric("mgnet_longer"))


setMethod("mgnet_longer", "mgnet", function(object,
                                            abundance=NULL, rel_abundance=NULL, norm_abundance=NULL,
                                            lineage=NULL, info_sample=NULL, info_taxa=NULL, community=NULL){
  

  # Checks
  if( isTRUE(abundance) & length(object@abundance) == 0 ) stop("abundance matrix missing.")
  if( isTRUE(rel_abundance) & length(object@rel_abundance) == 0 ) stop("rel_abundance matrix cannot be calculated, abundance or sample_sum missing. '?update_sample_sum' for more information.")
  if( isTRUE(norm_abundance) & length(object@norm_abundance) == 0 ) stop("norm_abundance matrix missing")
  if( isTRUE(lineage) & length(object@lineage) == 0) stop("lineage matrix missing")
  if( isTRUE(info_sample) & length(object@info_sample) == 0) stop("info_sample data.frame missing")
  if( isTRUE(info_taxa) & length(object@info_taxa) == 0) stop("info_taxa data.frame missing")
  if( isTRUE(community) & length(object@community) == 0) stop("community slot missing")

  # Automatic Assignments
  if( is.null(abundance) & length(object@abundance) != 0) abundance <- TRUE else abundance <- FALSE
  if( is.null(rel_abundance) & length(object@rel_abundance) != 0 ) rel_abundance <- TRUE else rel_abundance <- FALSE
  if( is.null(norm_abundance) & length(object@norm_abundance) != 0) norm_abundance <- TRUE else norm_abundance <- FALSE
  if( is.null(lineage) & length(object@lineage) != 0) lineage <- TRUE else lineage <- FALSE
  if( is.null(info_sample) & length(object@info_sample) != 0) info_sample <- TRUE else info_sample <- FALSE
  if( is.null(info_taxa) & length(object@info_taxa) != 0) info_taxa <- TRUE else info_taxa <- FALSE
  if( is.null(community) & length(object@community) != 0) community <- TRUE else community <- FALSE

  # reshape to longer format
  #----------------------------------------------------------------------------#
  
  # empty
  if( nsample(object)==0 & ntaxa(object)==0 ) return(tibble::tibble())

  # sample>0 and taxa>0
  if( nsample(object)!=0 & ntaxa(object)!=0) {
    
    long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                     taxa_id = taxa_id(object))
    
    if(abundance) {
      long_mgnet <- long_mgnet %>% left_join(
        abundance(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abundance"),
        by=c("sample_id","taxa_id"))
    }
    
    if(rel_abundance) {
      long_mgnet <- long_mgnet %>% left_join(
        rel_abundance(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "rel_abundance"),
        by=c("sample_id","taxa_id"))
    }
    
    if(norm_abundance) {
      long_mgnet <- long_mgnet %>% left_join(
        norm_abundance(object, .fmt="tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "norm_abundance"),
        by=c("sample_id","taxa_id"))
    }
    
    if(info_sample) {
      long_mgnet <- long_mgnet %>% left_join(
        info_sample(object, .fmt="tbl"), by="sample_id"
      )}
    
    if(lineage) {
      long_mgnet <- long_mgnet %>% left_join(
        lineage(object, .fmt="tbl"), by="taxa_id"
      )}
    
    if(info_taxa) {
      long_mgnet <- long_mgnet %>% left_join(
        info_taxa(object, .fmt="tbl"), by="taxa_id"
      )}
    
    if(community) {
      long_mgnet <- long_mgnet %>% left_join(
        community_members(object) %>% enframe( name = "taxa_id", value = "community"), 
        by = "taxa_id"
      )}
    
    return(long_mgnet)
    
  }
  
  # sample==0 taxa>0
  if( nsample(object)==0 & ntaxa(object)!=0) {
    
    long_mgnet <- tibble(taxa_id = taxa_id(object))
    
    if(lineage) {
      long_mgnet <- long_mgnet %>% left_join(
        lineage(object, .fmt="tbl"), by="taxa_id"
      )}
    
    if(info_taxa) {
      long_mgnet <- long_mgnet %>% left_join(
        info_taxa(object, .fmt="tbl"), by="taxa_id"
      )}
    
    return(long_mgnet)
    
  }
  
  # sample>0 taxa==0
  if( nsample(object)!=0 & ntaxa(object)==0) {
    
    long_mgnet <- tibble(taxa_id = sample_id(object))

    long_mgnet <- info_sample(object, "tbl")
    
    return(long_mgnet)
    
  }
  
})

setMethod("mgnet_longer", "mgnetList", function(object, source_mgnet_to,
                                                abundance=NULL, rel_abundance=NULL, norm_abundance=NULL,
                                                lineage=NULL, info_sample=NULL, info_taxa=NULL, community=NULL) {
  
  # Ensure source_mgnet_to is provided
  if(missing(source_mgnet_to) || source_mgnet_to == "") {
    stop("source_mgnet_to must be provided for mgnetList objects.")
  }
  
  # Initialize a list to store the longer format dataframes from each mgnet object
  longer_dfs <- list()
  
  # Make sure you're using the correct variable, `object`, which is the actual parameter name
  names_list <- names(object)  # Use 'object', not 'mgnetList'
  
  for (i in seq_along(object)) {  # Iterate over 'object', which is your actual mgnetList
    # Convert the current mgnet object to longer format
    longer_df <- mgnet_longer(object[[i]], 
                              abundance = abundance, rel_abundance = rel_abundance, norm_abundance = norm_abundance,
                              lineage = lineage, info_sample = info_sample, info_taxa = info_taxa, community = community)
    
    # Add a column to identify the source mgnet object, if not empty
    if (nrow(longer_df) > 0) {
      longer_df <- dplyr::mutate(longer_df, !!source_mgnet_to := names_list[i], .before = 1)  
      #longer_df <- dplyr::mutate(longer_df, source_mgnet = names_list[i], .before = 1)
    }
    
    # Store the result in the list
    longer_dfs[[i]] <- longer_df
  }
  
  # Combine all the dataframes into one, handling missing columns by filling with NA
  combined_df <- dplyr::bind_rows(longer_dfs)
  
  return(combined_df)
})