
#' # SET LINEAGE COLOR
#' #------------------------------------------------------------------------------#
#' #' Set Color Mapping to Top N Taxa in an `mgnet` Object
#' #'
#' #' This function assigns distinctive color codes to the top `n` taxa based on an aggregate metric computed from a specified data field
#' #' in an `mgnet` object. For `mgnetList` objects, it applies this logic individually to each contained `mgnet` object and
#' #' aggregates a unique set of top `n` taxa across all objects to ensure distinctive coloring. Taxa not in the top set across the list
#' #' are assigned a default or specified `extraColor`.
#' #'
#' #' @param object An `mgnet` or `mgnetList` object.
#' #' @param n Integer; the number of top taxa to be distinctly colored.
#' #' @param field Character; specifies the data field (`"abundance"`, `"rel_abundance"`, or `"norm_abundance"`) used for ranking taxa.
#' #' @param rank Optional character string specifying the taxonomic rank for aggregation and coloring.
#' #' @param order_fun Function to compute rankings, defaulting to `sum`. Other possibilities include `mean`, `median` or costum functions.
#' #' @param decreasing Logical indicating if sorting should be in decreasing order. When `TRUE` (default),
#' #' samples with the highest metric values are considered top.
#' #' @param extraColor String; hexadecimal color code for taxa outside the top set. This color is applied to taxa that do not make
#' #'        into the top `n` taxa across all `mgnet` objects in a `mgnetList` or within the single `mgnet` object. Default is "#FFFFFF" (white).
#' #' @param alpha Numeric; opacity level for colors, ranging from 0 (transparent) to 1 (opaque), with a default of 1.
#' #' @param colorspace Specifies the color space for generating the palette, with default 'pretty'. Options include 'pretty_dark', 'rainbow', 'pastels', or a custom list.
#' #' @param color_to Character; the name of the new column in `info_taxa` where color codes will be stored.
#' #' @return The updated `mgnet` or `mgnetList` object with `info_taxa` including color mappings.
#' #' @export
#' #' @importFrom dplyr mutate %>%
#' #' @name set_lineage_color
#' #' @aliases set_lineage_color,mgnet-method set_lineage_color,mgnetList-method
#' #' @seealso \link[mgnet]{colormap_categories}, \link[mgnet]{top_taxa}
#' setGeneric("set_lineage_color", function(object, n, 
#'                                          field, rank = NULL, order_fun = sum, decreasing = TRUE,
#'                                          extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
#'                                          color_to) standardGeneric("set_lineage_color"))
#' 
#' setMethod("set_lineage_color", signature = "mgnet", function(object, n, 
#'                                                              field, rank = NULL, order_fun = sum, decreasing = TRUE,
#'                                                              extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
#'                                                              color_to) {
#'   # Determine the top n taxa based on the specified criteria
#'   top_rank <- top_taxa(object, field, rank, order_fun, decreasing = decreasing)[1:n]
#'   all_rank <- unique(taxa_name(object,rank))
#'   
#'   categories <- c(top_rank, all_rank)
#'   categories <- categories[!duplicated(categories)]
#'   
#'   # Generate colors for the top taxa
#'   palette <- colormap_categories(categories = categories, distinctColors = n, extraColor = extraColor, alpha = alpha, colorspace = colorspace)
#'   
#'   object <- object %>%
#'     mutate_info_taxa(!!color_to := palette[taxa_name(object,rank)])
#'   
#'   validObject(object)
#'   return(object)
#' })
#' 
#' setMethod("set_lineage_color", signature = "mgnetList", function(object, n, 
#'                                                                  field, rank = NULL, order_fun = sum, decreasing = TRUE,
#'                                                                  extraColor = "#FFFFFF", alpha = 1, colorspace = "pretty",
#'                                                                  color_to) {
#'   # Determine the top n taxa based on the specified criteria
#'   top_rank <- top_taxa(object, field, rank, order_fun, decreasing = decreasing)[1:n]
#'   top_rank <- top_rank %>% unlist %>% unique
#'   n <- length(top_rank)
#'   
#'   if(n > 99) stop("the unique different taxa from all the different mgnet objects cannot be more than 99")
#'   all_rank <- unique(taxa_name(object,rank) %>% unlist %>% unique)
#'   
#'   categories <- c(top_rank, all_rank)
#'   categories <- categories[!duplicated(categories)]
#'   
#'   # Generate colors for the top taxa
#'   palette <- colormap_categories(categories = categories, distinctColors = n,
#'                                  extraColor = extraColor, alpha = alpha, colorspace = colorspace)
#'   
#'   object@mgnets <- sapply(object, \(x){
#'     x <- x %>%
#'       mutate_info_taxa(!!color_to := palette[taxa_name(x,rank)])
#'     return(x)
#'   })
#'   
#'   validObject(object)
#'   return(object)
#' })
#' 
#' 
#' # SET COMMUNITY COLOR
#' #------------------------------------------------------------------------------#
#' #' Set Color Mapping to Community in an `mgnet` Object or `mgnetList`
#' #'
#' #' This function assigns a visually distinct color palette to the communities detected in an `mgnet` object 
#' #' or across multiple `mgnet` objects within an `mgnetList`. It provides options to handle communities based on their sizes,
#' #' assigning unique colors to larger communities while using a common color for smaller ones, thus facilitating clearer
#' #' visual distinctions in subsequent analyses or visualizations.
#' #'
#' #' @param object An `mgnet` or `mgnetList` object.
#' #' @param sizes_threshold An integer specifying the minimum size a community must have to be assigned a unique color.
#' #'        Communities with sizes below this threshold will be colored using `smaller_color`. 
#' #' @param smaller_color A hexadecimal color code representing the color to be used for smaller communities, i.e., 
#' #'        those communities whose sizes do not meet the `sizes_threshold`. Defaults to light gray ("#D3D3D3").
#' #' @param isolated_color A hexadecimal color code used to represent isolated nodes, typically assigned to community '0' in mgnet.
#' #'        Defaults to white ("#FFFFFF"). See \link[mgnet]{as_mgnet_communities}.
#' #' @param alpha A numeric value between 0 and 1 that sets the transparency level of the colors, where 1 is fully opaque and
#' #'        0 is fully transparent.
#' #' @param colorspace Defines the scheme used to generate the color palette. Possible options are 'pretty', 'pretty_dark',
#' #'        'rainbow', 'pastels', or a custom list specifying the 'h' (hue), 's' (saturation), and 'l' (lightness) components.
#' #'        Each option provides a different aesthetic and can be chosen based on the visual requirements of the analysis.
#' #'        The function encapsulate qualpal to generate qualitative distinct color, see \link[qualpalr]{qualpal} for more details.
#' #' @param color_to The name of the new column within `info_taxa` where the assigned colors will be stored.
#' #' @param distinct_color Logical indicating whether distinct colors should be used for each `mgnet` object within an `mgnetList`
#' #'        when `object` is an `mgnetList`. If TRUE, each `mgnet` object's communities are assigned unique colors independently;
#' #'        if FALSE, a shared palette is used for all communities across the `mgnetList`, based on the maximum community number.
#' #'
#' #' @details
#' #' The function first assesses the community structure within the provided `mgnet` or `mgnetList` object. It applies a color
#' #' generation process that respects the specified `sizes_threshold`, assigning distinct colors to larger communities and a common
#' #' color to smaller ones. The `isolated_color` is used specifically for nodes that do not belong to any detected community or are
#' #' considered outliers. This function is particularly useful in network analyses where community detection results need to be
#' #' visualized to understand the structure and composition of the network.
#' #'
#' #' @return An updated `mgnet` or `mgnetList` object with the `info_taxa` slot modified to include a new column named as specified
#' #'         by `color_to`, which contains the color codes assigned based on community membership and size criteria.
#' #'
#' #' @export
#' #' @importFrom dplyr mutate %>%
#' #' @importFrom qualpalr qualpal
#' #' @importFrom grDevices adjustcolor
#' #' @name set_community_color
#' #' @aliases set_community_color,mgnet-method set_community_color,mgnetList-method
#' #' @seealso \link[mgnet]{colormap_communities} \link[mgnet]{cluster_signed}
#' setGeneric("set_community_color", function(object, 
#'                                            sizes_threshold, smaller_color = "#D3D3D3",
#'                                            isolated_color = "#FFFFFF", 
#'                                            alpha = 1, colorspace = "pretty",
#'                                            color_to, distinct_color = TRUE) standardGeneric("set_community_color"))
#' 
#' setMethod("set_community_color", signature = "mgnet", function(object, 
#'                                                                sizes_threshold, smaller_color = "#D3D3D3",
#'                                                                isolated_color = "#FFFFFF", 
#'                                                                alpha = 1, colorspace = "pretty",
#'                                                                color_to, distinct_color = TRUE){
#'   if(missing(sizes_threshold)){
#'     sizes_threshold <- 1
#'   } else {
#'     if(!is.numeric(sizes_threshold) || sizes_threshold != round(sizes_threshold) || sizes_threshold < 1) {
#'       stop("sizes_threshold must be a positive integer >= 1.")
#'     }
#'   }
#'   
#'   ncomm <- max(as.numeric(names(sizes(as_mgnet_communities(community(object))))))
#'   
#'   n <- igraph::sizes(as_igraph_communities(community(object)))
#'   n <- n[n>sizes_threshold]
#'   n <- max(as.numeric(names(n)))
#'   
#'   if(n > 99) {
#'     stop("There more than 99 distinct communities where i have to associate a distinct color. Unlucky, this function cannot generate more than 99 distinct color, try to increase the sizes_threshold to put in the same smaller communities under the color smaller_color")
#'   }
#'   
#'   if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
#'     stop("alpha must be in the range [0,1].")
#'   }
#'   
#'   if(!is.character(colorspace) && !is.list(colorspace)) {
#'     stop("colorspace must be a character or a list.")
#'   }
#'   
#'   if(is.character(colorspace)){
#'     if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
#'       stop("When colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
#'     }
#'   }
#'   
#'   if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
#'     stop("When colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpal.")
#'   }
#'   
#'   palette <- colormap_communities(community(object), 
#'                                   sizes_threshold = sizes_threshold, smaller_color = smaller_color,
#'                                   isolated_color = isolated_color, 
#'                                   colorspace = colorspace, alpha = alpha)
#'   
#'   object <- object %>%
#'     mutate_info_taxa(!!color_to := palette[community_members(object)])
#'   
#'   validObject(object)
#'   return(object)
#' })
#' 
#' setMethod("set_community_color", signature = "mgnetList", function(object, 
#'                                                                    sizes_threshold, smaller_color = "#D3D3D3",
#'                                                                    isolated_color = "#FFFFFF", 
#'                                                                    alpha = 1, colorspace = "pretty",
#'                                                                    color_to, distinct_color = TRUE){
#'   
#'   if(!missing(sizes_threshold)){
#'     if(!is.numeric(sizes_threshold) || sizes_threshold != round(sizes_threshold) || sizes_threshold <= 1) {
#'       stop("sizes_threshold must be a positive integer >= 1.")
#'     }
#'   }
#'   
#'   if(missing(sizes_threshold)){
#'     sizes_threshold <- 1
#'   }
#'   
#'   if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
#'     stop("alpha must be in the range [0,1].")
#'   }
#'   
#'   if(!is.character(colorspace) && !is.list(colorspace)) {
#'     stop("colorspace must be a character or a list.")
#'   }
#'   
#'   if(is.character(colorspace)){
#'     if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
#'       stop("When colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
#'     }
#'   }
#'   
#'   if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
#'     stop("When colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpal.")
#'   }
#'   
#'   if(distinct_color){
#'     
#'     ncommunity_sizes_threshold <- sapply(community(object), \(x){
#'       communities_sizes <- igraph::sizes(x)
#'       communities_sizes <- communities_sizes[communities_sizes >= sizes_threshold]
#'       return(max(as.numeric(names(communities_sizes))))
#'     })
#'     
#'     n_distinct_color <- sum(ncommunity_sizes_threshold)
#'     
#'     if(n_distinct_color > 99) {
#'       stop("There more than 99 distinct communities where i have to associate a distinct color. Unlucky, this function cannot generate more than 99 distinct color, try to increase the sizes_threshold to put in the same smaller communities under the color smaller_color")
#'     }
#'     
#'     colors <- qualpalr::qualpal(n = n_distinct_color, colorspace = colorspace)$hex
#'     colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
#'     
#'     for(i in 1:length(object)){
#'       
#'       mgnetObj <- object[[i]]
#'       
#'       start_color_idx <- ifelse(i==1, 1, sum(ncommunity_sizes_threshold[1:(i-1)])+1 )
#'       final_color_idx <- sum(ncommunity_sizes_threshold[1:i])
#'       colors_i <- colors[start_color_idx:final_color_idx]
#'       
#'       palette <- rep(NA_character_, length(community(mgnetObj)))
#'       names(palette) <- as.character(1:length(palette))
#'       palette[1:length(colors_i)] <- colors_i
#'       palette[is.na(palette)] <- smaller_color
#'       palette <- c("0" = isolated_color, palette)
#'       
#'       object@mgnets[[i]] <- mgnetObj %>%
#'         mutate_info_taxa(!!color_to := palette[community_members(mgnetObj)])
#'     }
#'     
#'   } else { # false distinct color
#'     
#'     ncommunity_sizes_threshold <- sapply(community(object), \(x){
#'       communities_sizes <- igraph::sizes(x)
#'       communities_sizes <- communities_sizes[communities_sizes >= sizes_threshold]
#'       return(max(as.numeric(names(communities_sizes))))
#'     })
#'     
#'     n_distinct_color <- max(ncommunity_sizes_threshold)
#'     
#'     if(n_distinct_color > 99) {
#'       stop("There more than 99 distinct communities where i have to associate a distinct color. Unlucky, this function cannot generate more than 99 distinct color, try to increase the sizes_threshold to put in the same smaller communities under the color smaller_color")
#'     }
#'     
#'     colors <- qualpalr::qualpal(n = n_distinct_color, colorspace = colorspace)$hex
#'     colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
#'     
#'     for(i in 1:length(object)){
#'       
#'       mgnetObj <- object[[i]]
#'       
#'       community_number <- ncommunity_sizes_threshold[i]
#'       colors_i <- colors[1:final_color_idx]
#'       
#'       palette <- rep(NA_character_, length(community(mgnetObj)))
#'       names(palette) <- as.character(1:length(palette))
#'       palette[1:length(colors_i)] <- colors_i
#'       palette[is.na(palette)] <- smaller_color
#'       palette <- c("0" = isolated_color, palette)
#'       
#'       object@mgnets[[i]] <- mgnetObj %>%
#'         mutate_info_taxa(!!color_to := palette[community_members(mgnetObj)])
#'     }
#'   }
#'   
#'   validObject(object)
#'   return(object)
#' })