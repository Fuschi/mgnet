#' Generate Colors Based on Factor Levels
#'
#' Simplified method to assign colors based on the levels of a pre-existing factor variable in the taxa metadata of `mgnet` or `mgnetList` objects.
#' 
#' @description
#' This function generates a color palette based on the levels of a specified factor variable in the taxa metadata. The colors are then stored as a new column.
#'
#' @param object An `mgnet` or `mgnetList` object containing taxa data with a factor variable.
#' @param var The name of the column (as a string) in the taxa metadata representing the factor variable.
#' @param color_to The name of the column (as a string) where the generated colors will be stored.
#' @param colorspace Defines the color space for color generation. It can be one of the predefined 
#'        character options ('pretty', 'pretty_dark', 'rainbow', 'pastels') or a list specifying 
#'        the 'h', 's', and 'l' elements for custom color spaces. For a detailed explanation of 
#'        these options, see the `qualpalr` package documentation or enter `??qualpalr` in the R console.
#' @param alpha Numeric; the transparency level of the colors, where 1 is opaque and 0 is fully transparent.
#'
#' @details
#' The function relies on the pre-existing factor levels in the specified column to generate distinct colors. 
#' This makes it easy to apply a consistent color scheme to ordered or lumped categories.
#'
#' @return The modified `mgnet` or `mgnetList` object with the new color column added.
#'
#' @export
#' @aliases colorize_taxa,mgnet-method colorize_taxa,mgnetList-method
setGeneric("colorize_taxa", function(object, var, color_to = NULL,
                                              colorspace = "pretty", alpha = 1) {
  standardGeneric("colorize_taxa")
})

setMethod("colorize_taxa", "mgnet", function(object, var, color_to = NULL,
                                                      colorspace = "pretty", alpha = 1) {
  
  if(miss_taxa(object)) stop("Error: no taxa available.")
  if(miss_metataxa(object)) stop("Error: no metadata on taxa available.")
  
  # Capture the variable as a quosure
  var_name <- rlang::as_string(rlang::quo_get_expr(rlang::enquo(var)))
  
  if (is.null(color_to)) {
    color_to <- paste0("color_", var_name)
  } else if (!is.character(color_to)) {
    stop("Error: 'color_to' must be a valid character string.")
  }
  
  var <- gather_taxa(object)[[var_name]]
  # Check if the column is a factor
  if (!is.factor(var)) {
    stop(paste("Error: The variable", var_name, "is not a factor in the taxa metadata."))
  }
  
  # Check if the column is a factor
  if (!is.factor(var)) {
    stop(paste("Error: The variable", var_name, "is not a factor in the taxa metadata."))
  }
  
  # Extract the levels from the factor variable
  factor_levels <- levels(var)
  
  # Generate color palette
  color_palette <- qualpalr::qualpal(length(factor_levels), colorspace = colorspace)
  color_palette <- rownames(color_palette$RGB)
  color_palette <- grDevices::adjustcolor(color_palette, alpha.f = alpha)
  names(color_palette) <- factor_levels
  # Add the color column to the taxa metadata
  taxa(object) <- gather_taxa(object) %>%
    dplyr::mutate(!!color_to := color_palette[as.character(.[[var_name]])]) 
  
  return(object)
})

setMethod("colorize_taxa", "mgnetList", function(object, var, color_to = NULL,
                                             colorspace = "pretty", alpha = 1) {
  
  if(miss_taxa(object, "any")) stop("Error: no taxa available in at least one mgnet element.")
  if(miss_metataxa(object, "any")) stop("Error: no metadata on taxa available.")
  
  # Capture the variable as a quosure
  var_name <- rlang::as_string(rlang::quo_get_expr(rlang::enquo(var)))
  
  if (is.null(color_to)) {
    color_to <- paste0("color_", var_name)
  } else if (!is.character(color_to)) {
    stop("Error: 'color_to' must be a valid character string.")
  }
  
  # Check if the factor variable exists
  if (!var_name %in% taxa_vars(object, "unique")) {
    stop(paste("Error: The factor variable", var, "does not exist in the taxa metadata."))
  }
  
  var <- gather_taxa(object)[[var_name]]
  # Check if the column is a factor
  if (!is.factor(var)) {
    stop(paste("Error: The variable", var_name, "is not a factor in the taxa metadata."))
  }
  
  # Extract the levels from the factor variable
  factor_levels <- levels(var)
  
  # Generate color palette
  color_palette <- qualpalr::qualpal(length(factor_levels), colorspace = colorspace)
  color_palette <- rownames(color_palette$RGB)
  color_palette <- grDevices::adjustcolor(color_palette, alpha.f = alpha)
  names(color_palette) <- factor_levels
  
  colors <- color_palette[as.character(var)]
  # Add the color column to the taxa metadata
  taxa(object) <- gather_taxa(object) %>%
    dplyr::mutate(!!rlang::ensym(color_to) := colors) %>%
    split_arrange_merged_taxa(object)
  
  return(object)
})


# SET COMMUNITIES COLOR
#------------------------------------------------------------------------------#
#' Set Community Colors in `mgnet` Objects
#'
#' This function applies a color scheme to taxa based on community membership, which is determined by 
#' the sizes of communities identified in a metagenomic network. Each community can be assigned a unique 
#' color, while smaller communities and isolated nodes are assigned default colors to maintain visual 
#' clarity.
#'
#' @description
#' `set_community_color` assigns colors to taxa based on community memberships that are specified in the `comm_id` 
#' of the `mgnet` object. The function employs a color palette that is sensitive to community size, differentiating 
#' between larger and smaller communities.
#'
#' @param object An `mgnet` object that includes community data.
#' @param size An integer that specifies the minimum size a community must have to receive a unique color.
#'        Communities smaller than this size will be colored with `smaller_color` (default 1).
#' @param color_to The name of the column in the taxa data where the color values will be stored (default 'color_comm').
#' @param smaller_color A hexadecimal color code used for communities that do not exceed the `size` threshold,
#'        defaulting to a light gray ("#D3D3D3").
#' @param isolated_color A hexadecimal color code used for isolated nodes, often those not belonging to any community,
#'        defaulting to white ("#FFFFFF").
#' @param colorspace The qualitative color palette from which to generate the community colors. Options include
#'        'pretty', 'pretty_dark', 'rainbow', 'pastels', and others as supported by the `qualpalr` package. This
#'        parameter allows customization of the color scheme to fit different visualization needs.
#' @param alpha A numeric value between 0 and 1 indicating the opacity of the colors, where 1 is fully opaque and
#'        0 is completely transparent. This is useful for creating layered visual effects in plots.
#'
#' @details
#' The function modifies the provided `mgnet` object by adding or updating a column with community colors.
#' It relies on the `colormap_community` function to generate a suitable color palette based on community
#' sizes and the specified colorspace. The visualization can then be tailored to highlight community structures
#' within metagenomic networks effectively.
#'
#' @return The modified `mgnet` object with the added or updated community color data.
#'
#' @export
#' @seealso \link[colormap_community]{Generate Color Palette for Community Visualization}
#' @aliases set_community_color,mgnet-method set_community_color,mgnetList-method
setGeneric("set_community_color", function(object, size, color_to = "color_comm",
                                           smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                           colorspace = "pretty", alpha = 1) standardGeneric("set_community_color"))

setMethod("set_community_color", signature = "mgnet", function(object, size, color_to = "color_comm",
                                                               smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                                               colorspace = "pretty", alpha = 1){
  
  if(miss_taxa(object)) stop("Error: No taxa available.")
  if(miss_slot(object, "comm")) stop("Error: No `comm` available.")
  
  palette <- colormap_community(comm(object), sizes_threshold = size, 
                                smaller_color = smaller_color, isolated_color = isolated_color,
                                colorspace = colorspace, alpha = alpha)
  
  
  object <- mutate_taxa(object, !!color_to := palette[as.character(comm_id)])
  return(object)
})


setMethod("set_community_color", signature = "mgnetList", function(object, size, color_to = "color_comm",
                                                               smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                                               colorspace = "pretty", alpha = 1){
  
  if(miss_taxa(object, "any")) stop("Error: No taxa available in at least one of the mgnet objects.")
  if(miss_slot(object, "comm", "any")) stop("Error: No `comm` available in at least one of the mgnet objects.")
  
  palette <- colormap_communities(comm(object), sizes_threshold = size, 
                                  smaller_color = smaller_color, isolated_color = isolated_color,
                                  colorspace = colorspace, alpha = alpha)
  
  for(i in names(object)){
    object[[i]] <- mutate_taxa(object[[i]], !!color_to := palette[[i]][as.character(comm_id)])
  }
  
  return(object)
})

