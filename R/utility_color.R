# TAXONOMY COLORMAP
#------------------------------------------------------------------------------#
#' Create Color Palette for Taxonomy
#'
#' Generates a distinct color palette for taxonomic classifications. Each unique taxon is assigned 
#' a color up to a specified limit (`distinctColor`). Beyond this limit, `extraColor` is used for 
#' the remaining taxa. Duplicate taxa identifiers are filtered out to ensure uniqueness.
#'
#' @param taxaID A character vector of taxonomic identifiers.
#' @param distinctColor An integer specifying the maximum number of distinct colors to generate, 
#'        ranging from 1 to 99.
#' @param extraColor The color used for taxa identifiers exceeding the `distinctColor` limit, 
#'        defaults to light gray. Accepts hexadecimal color codes.
#' @param alpha A numeric value specifying the transparency level of the colors.
#' @param colorspace Defines the color space for color generation. It can be one of the predefined 
#'        character options ('pretty', 'pretty_dark', 'rainbow', 'pastels') or a list specifying 
#'        the 'h', 's', and 'l' elements for custom color spaces. For a detailed explanation of 
#'        these options, see the `qualpalr` package documentation or enter `??qualpal` in the R console.
#'
#' @return A named vector of hexadecimal color codes, with names corresponding to the taxa IDs.
#'
#' @export
#' @importFrom qualpalr qualpal
#' @importFrom grDevices adjustcolor
colormap_taxonomy <- function(taxaID, distinctColor = 20, extraColor = "#BFBFBF",
                              alpha = 1, colorspace = "pretty") {
  
  if(!is.character(taxaID)) {
    stop("taxaID must be a character vector.")
  }
  
  taxaID <- unique(taxaID) # Select unique taxa IDs to avoid duplicate colors.
  
  if(!is.numeric(distinctColor) || distinctColor != round(distinctColor) || distinctColor < 1 || distinctColor > 99) {
    stop("distinctColor must be an integer between 1 and 99.")
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
  
  colors <- qualpalr::qualpal(n = distinctColor, colorspace = colorspace)$hex
  colors <- adjustcolor(colors, alpha.f = alpha)
  
  if(length(taxaID) > distinctColor) {
    extraColors <- rep(extraColor, length(taxaID) - distinctColor)
    extraColors <- adjustcolor(extraColors, alpha.f = alpha)
    colors <- c(colors, extraColors)
  }
  
  names(colors) <- taxaID
  
  return(colors)
}


# TAXONOMY COMMUNITY
#------------------------------------------------------------------------------#
#' Create Color Palette for Communities
#'
#' Generates a color palette with a specified number of distinct colors for visualizing communities. 
#' The first color, customizable via the `isolated_color` parameter, is designated for isolated nodes, 
#' with subsequent colors generated for other communities. This function uses the `qualpalr` package 
#' to create visually distinct colors.
#'
#' @param n A positive integer indicating the number of distinct colors to generate, not exceeding 99.
#' @param alpha A numeric value specifying the transparency level of the colors, where 1 is opaque 
#'        and 0 is fully transparent.
#' @param colorspace Defines the color space for color generation. It can be one of the predefined 
#'        character options ('pretty', 'pretty_dark', 'rainbow', 'pastels') or a list specifying 
#'        the 'h', 's', and 'l' elements for custom color spaces. For a detailed explanation of 
#'        these options, consult the `qualpalr` package documentation or enter `??qualpal` in the R console.
#' @param isolated_color The color assigned to isolated nodes, represented by community '0'. 
#'        Defaults to white ("#FFFFFF"). Accepts hexadecimal color codes.
#'
#' @return A named vector of colors in hexadecimal format, with names from '0' to 'n', where '0' 
#'         represents isolated nodes. The color for isolated nodes can be customized with the 
#'         `isolated_color` parameter.
#'
#' @export
#' @importFrom qualpalr qualpal
#' @importFrom grDevices adjustcolor
colormap_communities <- function(n = 20, alpha = 1, colorspace = "pretty",
                                 isolated_color = "#FFFFFF") {
  
  if(!is.numeric(n) || n != round(n) || n <= 0 || n > 99) {
    stop("n must be a positive integer less than or equal to 99.", call. = FALSE)
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
  colors <- qualpalr::qualpal(n = n + 1, colorspace = colorspace)$hex
  # Assign white to isolated nodes
  colors[1] <- isolated_color
  # Adjust color transparency
  colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
  
  # Assign names
  names(colors) <- as.character(0:n)
  
  return(colors)
}
