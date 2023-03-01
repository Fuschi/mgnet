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
#' @param n positive integer indicated the number of distinct color, n<100.
#' @param colorspace (Optional) Permits to set the color ranges. For more details
#' see \code{\link{qualpal}}.
#' @param alpha (Optional) Level of trasparency of colors.
#' 
#' @importFrom qualpalr qualpal
#' @importFrom grDevices rgb adjustcolor
#' @export
colormap_communities <- function(n=20, alpha=1,
                                 colorspace=list(h=c(0,360),s=c(.25,1),l=c(.25,.75))){
  
  if(round(n)!=n | !is.numeric(n) | n<=0 | n>99) stop("n must be an integer positive number (not greater or equal to 100)")
  if(!is.numeric(alpha)) stop("alpha must be numeric")
  if(alpha<0 | alpha>1) stop("alpha must be in range [0,1]")
  if(!is.list(colorspace)) stop("colorspace must be a list")
  if(any(names(list)!=c("h","s","l"))) stop("colorspace must have as elements h,s,l")
  
  colormap <- qualpalr::qualpal(n=n, colorspace=colorspace)
  colormap <- rownames(colormap$RGB)
  colormap <- c(grDevices::rgb(1,1,1),colormap)
  colormap <- grDevices::adjustcolor(colormap,alpha.f=alpha)
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
#' @param colorspace (Optional) Permits to set the color ranges. For more details
#' see \code{\link{qualpal}}.
#' @param alpha (Optional) Level of trasparency of colors.
#' 
#' @importFrom grDevices adjustcolor
#' @importFrom qualpalr qualpal
#' @export
colormap_taxonomy <- function(taxaID, alpha=1,
                              colorspace=list(h=c(0,360),s=c(.25,1),l=c(.5,.9))){
  
  if(!is.character(taxaID)) stop("taxaID must be character")
  if(!is.null(dim(taxaID))) stop("taxaID can't be a matrix")
  if(!is.numeric(alpha)) stop("alpha must be numeric")
  if(length(unique(taxaID))>99) stop("cannot generate more than 99 different color")
  if(alpha<0 | alpha>1) stop("alpha must be in range [0,1]")
  if(!is.list(colorspace)) stop("colorspace must be a list")
  if(any(names(list)!=c("h","s","l"))) stop("colorspace must have as elements h,s,l")
  
  taxaID <- unique(taxaID)
  n <- length(taxaID)
  
  colormap <- qualpalr::qualpal(n=n, colorspace)
  colormap <- rownames(colormap$RGB)
  colormap <- grDevices::adjustcolor(colormap, alpha.f=alpha)
  names(colormap) <- taxaID
  
  return(colormap)
}