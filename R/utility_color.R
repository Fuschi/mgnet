# CATEGORICAL COLOR MAP
#------------------------------------------------------------------------------#
#' Create Color Palette for Categorical Data
#'
#' Generates a distinct color palette for categorical data classifications. Each unique category 
#' is assigned a color up to a specified limit (`distinctColors`). Beyond this limit, `extraColor` 
#' is used for the remaining categories. Duplicate identifiers are filtered out to ensure uniqueness.
#'
#' @param categories A character vector of identifiers for the categories.
#' @param distinctColors An integer specifying the maximum number of distinct colors to generate, 
#'        ranging from 1 to 99. Default is the length of categories.
#' @param extraColor The color used for identifiers exceeding the `distinctColors` limit, 
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
colormap_categories <- function(categories, distinctColors, extraColor = "#FFFFFF",
                                 alpha = 1, colorspace = "pretty") {
  
  if(!is.character(categories)) {
    stop("categories must be a character vector.")
  }
  
  categories <- unique(categories) # Ensure category uniqueness.
  
  if(missing(distinctColors)){
    distinctColors <- length(categories)
  } 
  
  if(!is.numeric(distinctColors) || distinctColors != round(distinctColors) || distinctColors < 1 || distinctColors > 99){
    stop("distinctColors must be an integer between 1 and 99.")
  }
  
  if(distinctColors > min(99, length(categories))){
    stop("distinctColors cannot exceed the current number of categories: ", length(categories))
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
  colors <- qualpalr::qualpal(n = distinctColors, colorspace = colorspace)$hex
  colors <- adjustcolor(colors, alpha.f = alpha)
  
  # Add extra colors if necessary
  if(length(categories) > distinctColors) {
    extraColors <- rep(extraColor, length(categories) - distinctColors)
    extraColors <- adjustcolor(extraColors, alpha.f = alpha)
    colors <- c(colors, extraColors)
  }
  
  names(colors) <- categories
  
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
#' @return A named vector of colors in hexadecimal format. The vector's names correspond to community IDs,
#'         with '0' typically representing isolated nodes. Each community exceeding the size threshold is assigned
#'         a unique color from the generated palette, while others receive the `smaller_color`.
#'
#' @export
#' @importFrom qualpalr qualpal
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
  
  ncomm <- max(as.numeric(names(sizes(as_mgnet_communities(community)))))
  
  n <- igraph::sizes(as_igraph_communities(community))
  n <- n[n>sizes_threshold]
  n <- max(as.numeric(names(n)))
  
  if(n > 99) {
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
  
  # Generate colors using qualpalr
  palette_unique <- qualpalr::qualpal(n = n, colorspace = colorspace)$hex
  # Create palette with smaller community
  palette <- rep(NA_character_, length = ncomm)
  palette[1:length(palette_unique)] <- palette_unique
  palette[is.na(palette)] <- smaller_color
  # Assign white to isolated nodes
  palette <- c(isolated_color, palette)
  # Adjust color transparency
  palette <- grDevices::adjustcolor(palette, alpha.f = alpha)
  
  # Assign names
  names(palette) <- as.character(0:ncomm)
  
  return(palette)
}
