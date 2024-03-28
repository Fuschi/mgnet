# MGNET_LONGER
#------------------------------------------------------------------------------#
#' Convert mgnetList Object to Longer Format
#'
#' @description
#' Converts data from each `mgnet` object within an `mgnetList` into a longer format,
#' suitable for easy analysis and casting, by combining abundance, sample metadata,
#' taxonomic information, and other relevant data. A new column, `source_mgnet`, is added
#' to identify the source `mgnet` object for each row in the combined data.
#'
#' @usage mgnet_longer(object, abundance=NULL, relative=NULL, log_abundance=NULL,
#'        lineage=NULL, info_sample=NULL, info_taxa=NULL)
#'
#' @param object An `mgnetList` object.
#' @param abundance Logical; if TRUE, includes abundance data.
#' @param relative Logical; if TRUE, includes relative abundance data.
#' @param log_abundance Logical; if TRUE, includes log-transformed abundance data.
#' @param lineage Logical; if TRUE, includes lineage information.
#' @param info_sample Logical; if TRUE, includes sample metadata.
#' @param info_taxa Logical; if TRUE, includes taxonomic information.
#'
#' @return A `tibble` in a long format containing the combined data from each `mgnet` object
#' in the `mgnetList`, with an additional column `source_mgnet` indicating the source object.
#'
#' @export
#' @aliases mgnet_longer,mgnet-method mgnet_longer,mgnetList-method
#'
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer expand_grid
#' @importFrom dplyr select right_join bind_rows mutate
#' @importFrom rlang enquos
setGeneric("mgnet_longer", function(object,
                                    abundance=NULL, relative=NULL, log_abundance=NULL,
                                    lineage=NULL, info_sample=NULL, info_taxa=NULL) standardGeneric("mgnet_longer"))

setMethod("mgnet_longer", "mgnet", function(object,
                                            abundance=NULL, relative=NULL, log_abundance=NULL,
                                            lineage=NULL, info_sample=NULL, info_taxa=NULL){
  

  # Checks
  if( isTRUE(abundance) & length(object@abundance) == 0 ) stop("abundance matrix missing.")
  if( isTRUE(relative) & !check_sample_sum(object) | length(object@abundance) == 0 ) stop("relative matrix cannot be calculated, abundance or sample_sum missing. '?update_sample_sum' for more information.")
  if( isTRUE(log_abundance) & length(object@log_abundance) == 0 ) stop("log_abundance matrix missing")
  if( isTRUE(lineage) & length(object@lineage) == 0) stop("lineage matrix missing")
  if( isTRUE(info_sample) & length(object@info_sample) == 0) stop("info_sample data.frame missing")
  if( isTRUE(info_taxa) & length(object@info_taxa) == 0) stop("info_taxa data.frame missing")

  # Automatic Assignments
  if( is.null(abundance) & length(object@abundance) != 0) abundance <- TRUE else abundance <- FALSE
  if( is.null(relative) & length(object@abundance) != 0 & check_sample_sum(object)) relative <- TRUE else relative <- FALSE
  if( is.null(log_abundance) & length(object@log_abundance) != 0) log_abundance <- TRUE else log_abundance <- FALSE
  if( is.null(lineage) & length(object@lineage) != 0) lineage <- TRUE else lineage <- FALSE
  if( is.null(info_sample) & length(object@info_sample) != 0) info_sample <- TRUE else info_sample <- FALSE
  if( is.null(info_taxa) & length(object@info_taxa) != 0) info_taxa <- TRUE else info_taxa <- FALSE

  # reshape to longer format
  if( nsample(object)==0 & ntaxa(object)==0 ) return(tibble::tibble())
  if( nsample(object)!=0 & ntaxa(object)==0 ) stop("no samples found")
  if( nsample(object)==0 & ntaxa(object)!=0 ) stop("no taxa found")

  long_mgnet <- tidyr::expand_grid(sample_id = sample_id(object),
                                     taxa_id = taxa_id(object))

  if(abundance) {
    long_mgnet <- long_mgnet %>% left_join(
      abundance(object, .fmt="tbl") %>%
        tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abundance"),
      by=c("sample_id","taxa_id"))
  }

  if(relative) {
    long_mgnet <- long_mgnet %>% left_join(
      relative(object, .fmt="tbl") %>%
        tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "relative"),
      by=c("sample_id","taxa_id"))
  }

  if(log_abundance) {
    long_mgnet <- long_mgnet %>% left_join(
      log_abundance(object, .fmt="tbl") %>%
        tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "log_abundance"),
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

  return(long_mgnet)
})

setMethod("mgnet_longer", "mgnetList", function(object,
                                                abundance=NULL, relative=NULL, log_abundance=NULL,
                                                lineage=NULL, info_sample=NULL, info_taxa=NULL) {
  
  # Initialize a list to store the longer format dataframes from each mgnet object
  longer_dfs <- list()
  
  # Make sure you're using the correct variable, `object`, which is the actual parameter name
  names_list <- names(object)  # Use 'object', not 'mgnetList'
  
  for (i in seq_along(object)) {  # Iterate over 'object', which is your actual mgnetList
    # Convert the current mgnet object to longer format
    longer_df <- mgnet_longer(object[[i]], 
                              abundance = abundance, relative = relative, log_abundance = log_abundance,
                              lineage = lineage, info_sample = info_sample, info_taxa = info_taxa)
    
    # Add a column to identify the source mgnet object, if not empty
    if (nrow(longer_df) > 0) {
      longer_df <- dplyr::mutate(longer_df, source_mgnet = names_list[i], .before = 1)
    }
    
    # Store the result in the list
    longer_dfs[[i]] <- longer_df
  }
  
  # Combine all the dataframes into one, handling missing columns by filling with NA
  combined_df <- dplyr::bind_rows(longer_dfs)
  
  return(combined_df)
})