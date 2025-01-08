# PLOT MGNET
#------------------------------------------------------------------------------#
#' Plot mgnet Network
#'
#' Visualizes an `mgnet` network object, emphasizing the signed properties of edges by 
#' color-coding them and adjusting vertex sizes and edge widths for enhanced readability.
#' Positive edges are colored differently from negative edges, and vertex sizes can be dynamically
#' scaled based on their properties.
#'
#' @description
#' `plot.mgnet` focuses on the visualization of signed networks. It highlights the positive
#' connections while still allowing adjustments to graphical representations of both vertices and
#' edges. The function is highly customizable, supporting various igraph layout options and 
#' additional aesthetic parameters.
#'
#' Vertex sizes are dynamically calculated using the formula:
#' \deqn{v=\text{\footnotesize{sumConst}}+\text{\footnotesize{multConst}}\left(\text{\footnotesize{max}}(v_0)\left(\frac{v_0}{\text{\footnotesize{max}}(v_0)}\right)^{\text{expFactor}}\right)}
#' where \eqn{v_0} is the initial vertex size. This formula allows for the enhancement of 
#' size differences between vertices, making the network's structure more apparent.
#'
#' @param x An `mgnet` object containing the network to be plotted.
#' @param layout An igraph layout function or coordinates for vertex positions. 
#'        If missing, `layout_signed` is applied, focusing on the structure formed by positive edges.
#' @param ... Additional arguments passed to `igraph::plot.igraph`, allowing further 
#'        customization of the plot. This includes all parameters supported by `plot.igraph`, 
#'        such as `vertex.label`, `edge.arrow.size`, among others.
#' @param expFactor A numeric value that adjusts the relative sizes of vertices. 
#'        Higher values increase the size difference between smaller and larger vertices.
#' @param multConst A positive numeric value applied as a multiplicative factor to all vertex sizes.
#' @param sumConst A positive numeric value added to the size of all vertices, 
#'        allowing for baseline size adjustment.
#' @param thickness A logical value; when `TRUE`, scales edge widths based on their absolute weights.
#' @param alphaFactor A numeric value in \[0,1\], adjusting the transparency of edge colors.
#' @param widthFactor A positive numeric value that scales the width of all edges uniformly.
#' @param posCol The color used for positive edges, specified as an RGB value.
#' @param negCol The color used for negative edges, specified as an RGB value.
#' @param maxSize An optional positive integer specifying the maximum size for vertices. 
#'        If set, vertex sizes are scaled such that the largest vertex has this size.
#' @param maxWidth An optional positive integer specifying the maximum width for edges. 
#'        If set, edge widths are scaled such that the widest edge has this width.
#' @param force.positive A logical indicating whether to adjust vertex sizes to be strictly positive 
#'        by offsetting them by the minimum size value found in the network.
#'
#' @return Generates a plot of the specified `mgnet` network with the given aesthetic parameters.
#'
#' @importFrom igraph layout.fruchterman.reingold subgraph.edges E<- E V V<-
#' @importFrom grDevices rgb adjustcolor
#' @export
#' @rdname plot.mgnet
plot.mgnet <- function(x,layout,
                       ...,
                       expFactor=1, multConst=1, sumConst=0,
                       thickness=TRUE,
                       alphaFactor=.5, widthFactor=1,
                       posCol=rgb(0,0,1), negCol=rgb(1,0,0),
                       maxSize=NULL, maxWidth=NULL,
                       force.positive=TRUE
) {
  
  if(length(netw(x))==0) stop("missing network in mgnet")
  if(!is.numeric(expFactor)) stop("expFactor must be numeric")
  if(!is.numeric(multConst) | multConst<=0) stop("multFact must be a number greater than 0")
  if(!is.numeric(sumConst)) stop("sumConst must be a number")
  if(!is.logical(thickness)) stop("thickness must be logical")
  if(!is.numeric(alphaFactor) | alphaFactor<0 | alphaFactor>1) stop("alphaFactor must be a number in range [0,1]")
  if(!is.numeric(widthFactor) | widthFactor<=0) stop("widthFactor must be a positive number")
  if(!is.logical(force.positive)) stop("force.positive must be logical")
  
  if(missing(layout)){
    layout <- layout_signed(netw(x))
  }
  
  g <- netw(x)
  
  # Edges
  if(thickness){
    E(g)$width <- abs(E(g)$weight)
    E(g)$width <- widthFactor*E(g)$width
    if(!is.null(maxWidth)){
      E(g)$width <- (E(g)$width/max(E(g)$width))*maxWidth
    }
  } else {
    E(g)$width <- 1
  }
  E(g)$color <- ifelse(E(g)$weight>0,posCol,negCol)
  E(g)$color <- adjustcolor(E(g)$color,alpha.f=alphaFactor)
  
  # Nodes
  params <- list(...)
  if(!is.null(V(g)$size)){
    vertexSize <- V(g)$size
  }
  if("vertex.size"%in%names(params)){
    vertexSize <- params[["vertex.size"]]
  } else {
    vertexSize <- 15
  }
  
  if(force.positive & any(vertexSize<0)){
    vertexSize <- vertexSize - min(vertexSize)
  }
  
  vertexSize <- max(vertexSize)*(vertexSize/max(vertexSize))^(expFactor)
  vertexSize <- sumConst + multConst* vertexSize
  if(!is.null(maxSize)){
    vertexSize <- (vertexSize/max(vertexSize))*maxSize
  }
  
  plot(g, layout=layout, vertex.size=vertexSize, ...)
}