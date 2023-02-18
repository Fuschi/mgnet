################################################################################
################################################################################
# COMMUNITIES COLORMAP
################################################################################
################################################################################
#' Create color palette for communities
#' 
#' @description User wrapper for distinctColorPalette function of randomcoloR
#' package. Create n+1 distinct color to associate with communities ID. The 
#' first color labeled as 0 is always white for isolated nodes.
#' 
#' @param n positive integer indicated the number of distinct color.
#' 
#' @importFrom qualpalr qualpal
#' @importFrom grDevices rgb
#' @export
colormap_communities <- function(n=20){
  
  if(round(n)!=n | !is.numeric(n) | n<=0 | n>99) stop("n must be an integer positive number (not greater or equal to 100)")
  
  colormap <- qualpalr::qualpal(n=n)
  colormap <- rownames(colormap$RGB)
  colormap <- c(grDevices::rgb(1,1,1,.8),colormap)
  names(colormap) <- as.character(0:n)
  
  return(colormap)
}
################################################################################
################################################################################
# END COMMUNITIES COLORMAP
################################################################################
################################################################################




################################################################################
################################################################################
# TAXONOMY COLORMAP
################################################################################
################################################################################
#' Create color palette for taxonomy
#' 
#' @description User wrapper for distinctColorPalette function of randomcoloR
#' package. Create named vector with taxaID as name and colors as values. If are
#' present duplicated in the vector the function choice only unique values.
#' 
#' @param taxaID character vector with all taxonomic classification.
#' 
#' @importFrom qualpalr qualpal
#' @export
colormap_taxonomy <- function(taxaID){
  
  if(!is.character(taxaID)) stop("taxaID must be character")
  if(!is.null(dim(taxaID))) stop("taxaID can't be a matrix")
  
  taxaID <- unique(taxaID)
  n <- length(taxaID)
  
  colormap <- qualpalr::qualpal(n=n)
  colormap <- rownames(colormap$RGB)
  names(colormap) <- taxaID
  
  return(colormap)
}