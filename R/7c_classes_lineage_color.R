# Set Taxonomic Color Mapping in mgnet Objects
#------------------------------------------------------------------------------#
#' Set Color Mapping for a Specified Taxonomic Rank in an mgnet Object
#'
#' This function updates the `info_taxa` data frame of an `mgnet` object to include a new
#' column named `color`, which contains hexadecimal color codes assigned to each taxon
#' at the specified taxonomic rank. The colors are generated using the `colormap_categorical` function.
#'
#' @param mgnet_obj An `mgnet` object.
#' @param rank A character string specifying the taxonomic rank at which to assign colors.
#'        The rank must be a column name within the `lineage` data frame of the `mgnet` object.
#' @param colorspace The color space setting for generating the palette (defaults to "pretty").
#'        Other options include 'pretty_dark', 'rainbow', 'pastels', or a custom list specifying 
#'        color space parameters.
#' @return An updated `mgnet` object with color mapping added to the `info_taxa`.
#' @export
#' @importFrom dplyr mutate left_join
#' @importFrom methods validObject
setGeneric("set_lineage_color", function(object, rank, distinctColor, extraColor = "#FFFFFF", 
                                         alpha = 1, colorspace = "pretty", color_to) standardGeneric("set_lineage_color"))
setMethod("set_lineage_color", signature = "mgnet", function(object, rank, distinctColor, extraColor = "#FFFFFF", 
                                                             alpha = 1, colorspace = "pretty", color_to) {
  if(!rank %in% colnames(object@lineage)) {
    stop("The specified rank is not available in the lineage data.")
  }
  
  # Extract the specific rank data
  rankData <- taxa_name(object,rank)
  
  # Ensure distinctColor does not exceed the number of unique categories
  if(missing(distinctColor)) {
    distinctColor <- length(unique(rankData))
  }
  
  # Generate colors using the colormap_categories function
  palette <- colormap_categories(categories = rankData, distinctColor = distinctColor, extraColor = extraColor, alpha = alpha, colorspace = colorspace)
  
  # Add or update the color information in the info_taxa
  if(length(info_taxa(object)) == 0){
    object@info_taxa <- data.frame(color_to = palette[rankData], row.names = taxa_id(object))
  } else {
    object@info_taxa <- dplyr::mutate(object@info_taxa, !!color_to := palette[rankData])
  }
  
  # Return the updated object
  return(object)
})
