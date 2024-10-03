# CATEGORICAL COLOR MAP
#------------------------------------------------------------------------------#
#' Create Color Palette for Categorical Data
#'
#' Generates a distinct color palette for categorical data classifications. Each unique category 
#' is assigned a color up to a specified limit (`distinct_colors`). Beyond this limit, `extraColor` 
#' is used for the remaining categories. Duplicate identifiers are filtered out to ensure uniqueness.
#'
#' @param categories A character vector of identifiers for the categories.
#' @param distinct_colors An integer specifying the maximum number of distinct colors to generate, 
#'        ranging from 1 to 99. Default is the length of categories.
#' @param extraColor The color used for identifiers exceeding the `distinct_colors` limit, 
#'        defaults to white. Accepts hexadecimal color codes.
#' @param alpha A numeric value specifying the transparency level of the colors.
#' @param colorspace Defines the color space for color generation. It can be one of the predefined 
#'        character options ('pretty', 'pretty_dark', 'rainbow', 'pastels') or a list specifying 
#'        the 'h', 's', and 'l' elements for custom color spaces. For a detailed explanation of 
#'        these options, see the `qualpalr` package documentation or enter `??qualpalr` in the R console.
#'
#' @return A named vector of hexadecimal color codes, with names corresponding to the category IDs.
#'
#' @export
#' @importFrom qualpalr qualpal
#' @importFrom grDevices adjustcolor
colormap_categories <- function(categories, distinct_colors, extraColor = "#FFFFFF",
                                 alpha = 1, colorspace = "pretty") {
  
  if(!is.character(categories) & !is.logical(categories)) {
    stop("categories must be a character or logical vector.")
  }
  
  uniq_categories <- unique(categories) # Ensure category uniqueness.
  
  if(missing(distinct_colors)){
    distinct_colors <- length(uniq_categories)
  } 
  
  if(!is.numeric(distinct_colors) || distinct_colors != round(distinct_colors) || distinct_colors < 1 || distinct_colors > 99){
    stop("distinct_colors must be an integer between 1 and 99.")
  }
  
  if(distinct_colors > min(99, length(uniq_categories))){
    stop("distinct_colors cannot exceed the current number of categories: ", length(uniq_categories))
  }
  
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be in the range [0,1].")
  }
  
  if(!is.character(colorspace) && !is.list(colorspace)) {
    stop("colorspace must be a character or a list.")
  }
  
  if(is.character(colorspace)){
    if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
      stop("When colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
    }
  }
  
  if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
    stop("When colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpalr.")
  }
  
  # Generate colors
  colors <- qualpalr::qualpal(n = distinct_colors, colorspace = colorspace)$hex
  colors <- adjustcolor(colors, alpha.f = alpha)
  
  # Add extra colors if necessary
  if(length(uniq_categories) > distinct_colors) {
    extraColors <- rep(extraColor, length(uniq_categories) - distinct_colors)
    extraColors <- adjustcolor(extraColors, alpha.f = alpha)
    colors <- c(colors, extraColors)
  }
  
  names(colors) <- uniq_categories
  
  return(colors)
}


# COLORMAP COMMUNITY
#------------------------------------------------------------------------------#
#' Create Color Palette for Community Visualization
#'
#' Generates a visually distinct color palette for representing community memberships in network
#' visualizations. The palette is tailored to emphasize larger communities by providing unique colors 
#' for those that exceed a size threshold. Smaller communities and isolated nodes receive a common or 
#' default color to maintain visual clarity. This function uses the `qualpalr` package for color generation,
#' supporting up to 99 distinct communities due to practical limits on visual distinction.
#'
#' @param community An object representing community data, typically derived from a community detection
#'        method such as those available in the `igraph` package.
#' @param sizes_threshold An integer specifying the minimum size a community must exceed to receive a
#'        unique color. If it missing the function generate color for all the communities also with size 2.
#'        Communities not meeting this size threshold will be assigned the `smaller_color`.
#' @param smaller_color A hexadecimal color code used for communities that do not exceed the size threshold,
#'        defaulting to a light gray ("#D3D3D3").
#' @param isolated_color A hexadecimal color code used for isolated nodes, often represented as community '0',
#'        defaulting to white ("#FFFFFF"). 
#' @param colorspace Specifies the color space from which to generate the palette. Available options include
#'        'pretty', 'pretty_dark', 'rainbow', 'pastels', or a custom list detailing the 'h' (hue), 's' (saturation),
#'        and 'l' (lightness) components. For detailed settings, refer to the `qualpalr` documentation.
#' @param alpha A numeric value between 0 and 1 that specifies the opacity of the colors, where 1 is fully opaque
#'        and 0 is fully transparent.
#'
#'
#' @export
#' @importFrom qualpalr qualpal
#' @importFrom igraph sizes
#' @importFrom grDevices adjustcolor
colormap_communities <- function(community, 
                                 sizes_threshold, smaller_color = "#D3D3D3",
                                 isolated_color = "#FFFFFF",
                                 colorspace = "pretty", alpha = 1) {
  
  if (!inherits(community, "communities")) {
    stop("The 'community' parameter must be a community object from igraph.")
  }
  
  
  if(missing(sizes_threshold)){
    sizes_threshold <- 1
  } else {
    if(!is.numeric(sizes_threshold) || sizes_threshold != round(sizes_threshold) || sizes_threshold < 1) {
      stop("sizes_threshold must be a positive integer >= 1.")
    }
  }
  
  # Get vector of community sizes
  vec_sizes <- igraph::sizes(community)
  
  # Total number of communities
  n_tot_comm <- max(as.numeric(names(vec_sizes[vec_sizes > 1])))
  
  # Number of community without the smaller ones given by sizes_threshold
  n_thresh_comm <- max(as.numeric(names(vec_sizes[vec_sizes>sizes_threshold])))

  if(n_thresh_comm > 99) {
    stop("There more than 99 distinct communities where i have to associate a distinct color. Unlucky, this function cannot generate more than 99 distinct color, try to increase the sizes_threshold to put in the same smaller communities under the color smaller_color")
  }
  
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be in the range [0,1].")
  }
  
  if(!is.character(colorspace) && !is.list(colorspace)) {
    stop("colorspace must be a character or a list.")
  }
  
  if(is.character(colorspace)){
    if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
      stop("When colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
    }
  }
  
  if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
    stop("When colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpal.")
  }
  
  # Generate colors using qualpalr for larger communities > sizes_threshold
  palette_only_large <- qualpalr::qualpal(n = n_thresh_comm, colorspace = colorspace)$hex
  # Create palette with smaller community
  palette_with_smaller <- rep(NA_character_, length = n_tot_comm)
  palette_with_smaller[1:length(palette_only_large)] <- palette_only_large
  palette_with_smaller[is.na(palette_with_smaller)] <- smaller_color
  # Assign white to isolated nodes
  palette_with_isolated <- rep(NA_character_, length = length(igraph::sizes(community)))
  palette_with_isolated[1:length(palette_with_smaller)] <- palette_with_smaller
  palette_with_isolated[is.na(palette_with_isolated)] <- isolated_color
  # Adjust color transparency
  palette_with_isolated <- grDevices::adjustcolor(palette_with_isolated, alpha.f = alpha)
  
  # Assign names
  names(palette_with_isolated) <- as.character(names(igraph::sizes(community)))
  
  return(palette_with_isolated)
}
