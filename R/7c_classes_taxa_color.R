# SET LINEAGE COLOR
#------------------------------------------------------------------------------#
#' Set Color Mapping to Top N Taxa in an `mgnet` Object
#'
#' This function assigns distinctive color codes to the top `n` taxa based on an aggregate metric computed from a specified data field
#' in an `mgnet` object. For `mgnetList` objects, it applies this logic individually to each contained `mgnet` object and
#' aggregates a unique set of top `n` taxa across all objects to ensure distinctive coloring. Taxa not in the top set across the list
#' are assigned a default or specified `extraColor`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param n Integer; the number of top taxa to be distinctly colored.
#' @param field Character; specifies the data field (`"abundance"`, `"rel_abundance"`, or `"norm_abundance"`) used for ranking taxa.
#' @param rank Optional character string specifying the taxonomic rank for aggregation and coloring.
#' @param order_fun Function to compute rankings, defaulting to `sum`. Other possibilities include `mean`, `median` or custum functions.
#' @param extraColor String; hexadecimal color code for taxa outside the top set. This color is applied to taxa that do not make
#'        into the top `n` taxa across all `mgnet` objects in a `mgnetList` or within the single `mgnet` object. Default is "#FFFFFF" (white).
#' @param alpha Numeric; opacity level for colors, ranging from 0 (transparent) to 1 (opaque), with a default of 1.
#' @param colorspace Specifies the color space for generating the palette, with default 'pretty'. Options include 'pretty_dark', 'rainbow', 'pastels', or a custom list.
#' @param color_to Character; the name of the new column in `info_taxa` where color codes will be stored.
#' @return The updated `mgnet` or `mgnetList` object with `info_taxa` including color mappings.
#' @export
#' @importFrom dplyr mutate %>%
#' @name set_lineage_color
#' @aliases set_lineage_color,mgnet-method set_lineage_color,mgnetList-method
#' @seealso \link[mgnet]{colormap_categories}, \link[mgnet]{top_n_taxa}
setGeneric("set_lineage_color", function(object, n, field, rank = NULL, order_fun = sum,
                                         extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                         color_to) standardGeneric("set_lineage_color"))

setMethod("set_lineage_color", signature = "mgnet", function(object, n, field, rank = NULL, order_fun = sum,
                                                             extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                                             color_to) {
  # Determine the top n taxa based on the specified criteria
  top_rank <- top_n_taxa(object, n, field, rank, order_fun, decreasing = TRUE)
  all_rank <- unique(taxa_name(object,rank))
  
  categories <- c(top_rank, all_rank)
  categories <- categories[!duplicated(categories)]
  
  # Generate colors for the top taxa
  palette <- colormap_categories(categories = categories, distinctColors = n, extraColor = extraColor, alpha = alpha, colorspace = colorspace)
  
  object <- object %>%
    mutate_info_taxa(!!color_to := palette[taxa_name(object,rank)])
  
  validObject(object)
  return(object)
})

setMethod("set_lineage_color", signature = "mgnetList", function(object, n, field, rank = NULL, order_fun = sum,
                                                             extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                                             color_to) {
  # Determine the top n taxa based on the specified criteria
  top_rank <- top_n_taxa(object, n, field, rank, order_fun, decreasing = TRUE)
  top_rank <- top_rank %>% unlist %>% unique
  n <- length(top_rank)
  
  if(n > 99) stop("the unique different taxa from all the different mgnet objects cannot be more than 99")
  all_rank <- unique(taxa_name(object,rank) %>% unlist %>% unique)
  
  categories <- c(top_rank, all_rank)
  categories <- categories[!duplicated(categories)]
  
  # Generate colors for the top taxa
  palette <- colormap_categories(categories = categories, distinctColors = n,
                                 extraColor = extraColor, alpha = alpha, colorspace = colorspace)
  
  object@mgnets <- sapply(object, \(x){
    x <- x %>%
      mutate_info_taxa(!!color_to := palette[taxa_name(x,rank)])
    return(x)
  })
  
  validObject(object)
  return(object)
})


# SET COMMUNITY COLOR
#------------------------------------------------------------------------------#
#' Set Color Mapping to Community in an `mgnet` Object
#'
#' @export
#' @importFrom dplyr mutate %>%
#' @name set_community_color
#' @aliases set_community_color,mgnet-method
#' @seealso \link[mgnet]{colormap_communities}
setGeneric("set_community_color", function(object, n,
                                           isolated_color = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                           color_to, distinct_color = TRUE) standardGeneric("set_community_color"))

setMethod("set_community_color", signature = "mgnet", function(object, n,
                                                             isolated_color = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                                             color_to, distinct_color = TRUE){
  
  if(missing(n)) n <- max(membership(community(object)))
  if(n > 99) stop("the unique different taxa from all the different mgnet objects cannot be more than 99")
  palette <- colormap_communities(n = n, alpha = alpha, colorspace = colorspace, isolated_color = isolated_color)
  
  object <- object %>%
    mutate_info_taxa(!!color_to := palette[as_mgnet_communities(community(object))])
  
  return(object)
})

setMethod("set_community_color", signature = "mgnetList", function(object, n,
                                                               isolated_color = "#FFFFFF", alpha = 1, colorspace = "pretty",
                                                               color_to, distinct_color = NULL){
  
  if(is.null(distinct_color)){
    
  }
  
  
  if(n > 99) stop("the unique different taxa from all the different mgnet objects cannot be more than 99")
  palette <- colormap_communities(n = n, alpha = alpha, colorspace = colorspace, isolated_color = isolated_color)
  
  object <- object %>%
    mutate_info_taxa(!!color_to := palette[as_mgnet_communities(community(object))])
  
  return(object)
})

