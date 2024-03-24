# MGNET EXCTRACTOR
#------------------------------------------------------------------------------#
#' Subset mgnet Object
#'
#' Subsets an `mgnet` object based on provided sample and taxa indices. 
#' This method allows for the selective extraction of abundance, log abundance data, 
#' sample and taxa metadata, and optionally, network and community data based on 
#' specified indices.
#'
#' @param x An `mgnet` object to be subsetted.
#' @param i Indices or logical vector indicating samples to be included in the subset. If missing, all samples are included.
#' @param j Indices or logical vector indicating taxa to be included in the subset. If missing, all taxa are included.
#'
#' @details The method enables the extraction of specific parts of the `mgnet` 
#' object based on sample and taxa indices. When subsetting by samples, the abundance, 
#' log abundance, and sample metadata are filtered accordingly. 
#' Subsetting by taxa filters the taxa-related data and, if present, updates the 
#' network and community slots to reflect the selected taxa. 
#' It's important to note that subsetting by samples remove the network and community slots, 
#' and a warning will be sended in such cases.
#'
#' @return A new `mgnet` object containing only the specified subset of samples and taxa. The structure and type of data within the `mgnet` object are preserved.
#' @seealso \code{\link[base]{Extract}} for base subsetting functions.
#' 
#' @export
setMethod(f="[", signature="mgnet",function(x,i,j){
  
  ifelse(length(x@abundance)!=0, abundance.new<-x@abundance[i,j,drop=F], abundance.new<-x@abundance)
  ifelse(length(x@info_sample)!=0, info_sample.new<-x@info_sample[i, ,drop=F], info_sample.new<-x@info_sample)
  ifelse(length(x@lineage)!=0, lineage.new<-x@lineage[j, ,drop=F], lineage.new<-x@lineage)
  ifelse(length(x@info_taxa)!=0, info_taxa.new<-x@info_taxa[j, ,drop=F], info_taxa.new<-x@info_taxa)
  ifelse(length(x@log_abundance)!=0, log_abundance.new<-x@log_abundance[i,j,drop=F], log_abundance.new<-x@log_abundance)
  
  if(length(x@network)==0){
    
    return(mgnet(abundance=abundance.new,
                 info_sample=info_sample.new,
                 lineage=lineage.new,
                 info_taxa=info_taxa.new,
                 log_abundance=log_abundance.new))
    
    } else if(length(x@network)!=0 & missing(i)){
      
      ifelse(length(x@network)!=0,
             network.new<-igraph::subgraph(x@network,j),
             network.new<-x@network)
      
      if(length(x@community)!=0){
        community.new <- x@community
        if(is.character(j)) j <- which(taxa_id(x)%in%j)
        community.new$membership <- x@community$membership[j]
        community.new$vcount <- length(community.new$membership)
        community.new$modularity <- NA
      } else {
        community.new <- x@community
      }
      
      return(mgnet(abundance = abundance.new,
                   info_sample = info_sample.new,
                   lineage = lineage.new,
                   info_taxa = info_taxa.new,
                   log_abundance = log_abundance.new,
                   network = network.new,
                   community = community.new))
      
  } else {
    
    warning("sample subsetting is not defined for network the resulting object will not have the slots network and community")
    return(mgnet(abundance=abundance.new,
                 info_sample=info_sample.new,
                 lineage=lineage.new,
                 info_taxa=info_taxa.new,
                 log_abundance=log_abundance.new
    ))
  }
})



#' Extract Subsets from mgnetList
#'
#' This method enables the extraction of sub-lists or individual `mgnet` objects from an `mgnetList`
#' container by specifying numeric indices or names. Use `[` to extract sub-lists or individual
#' objects and maintain the `mgnetList` structure.
#'
#' @param x An `mgnetList` object.
#' @param i Numeric indices or character vector of names specifying the `mgnet` objects to extract.
#'          If missing, the entire `mgnetList` is returned.
#' @param j Ignored in this method, included for compatibility with the generic method.
#' @param drop Ignored in this method, included for compatibility with the generic method.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Depending on the type and number of indices or names provided in `i`:
#'         - A single `mgnet` object if a single index or name is provided.
#'         - An `mgnetList` containing the selected `mgnet` objects if multiple indices or names are provided.
#'         - The full `mgnetList` if `i` is missing.
#'
#' @seealso \code{\link[base]{Extract}} for base R extraction operations.
#' @export
setMethod(f="[", signature="mgnetList", function(x, i, j, ...) {
  
  if(!missing(j)) stop("wrong dimension number")
  
  if(missing(i)) {
    
    return(x)
    
  } else {
    
    if(is.numeric(i)) {
      if(any(i > length(x@mgnets)) || i < 1) {
        stop("Index out of bounds.")
      }
      return(mgnetList(x@mgnets[i]))
      
    } else if(is.character(i)) {
      
      if(!all(i %in% names(x@mgnets))) {
        stop("No mgnet object with such a name in the list.")
      }
      return(mgnetList(x@mgnets[i]))
      
    } else {
      stop("Index must be either numeric or a character vector of names.")
    }
  }
})


#' Extract mgnet Objects from mgnetList
#'
#' This method extracts one or more `mgnet` objects from an `mgnetList` using either numeric indices
#' or names. When a single index or name is provided, a single `mgnet` object is extracted. When multiple
#' indices or names are provided, a simple list of `mgnet` objects is returned.
#'
#' Use `[[` when you need to work directly with the extracted `mgnet` object(s). This method allows
#' for flexibility in extracting and working with individual or multiple `mgnet` objects outside
#' of an `mgnetList` structure.
#'
#' @param x An `mgnetList` object.
#' @param i Numeric index, vector of numeric indices, name, or vector of names specifying
#'          the `mgnet` object(s) to extract. If missing, an error is thrown.
#' @param j Ignored in this method, included for compatibility with the generic method.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Depending on the input `i`:
#'         - A single `mgnet` object if `i` specifies a single index or name.
#'         - A simple list of `mgnet` objects if `i` specifies multiple indices or names.
#' @seealso \code{\link[base]{Extract}} for base R extraction operations.
#' @export
setMethod(f="[[", signature="mgnetList", function(x, i, j, ...) {
  
  if(!missing(j)) stop("wrong dimension number")
  
  if(missing(i)) {
    
    return(x@mgnets)
    
  } else {
    
    if(is.numeric(i)) {
      if(any(i > length(x@mgnets)) || i < 1) {
        stop("Index out of bounds.")
      }
      return(x@mgnets[i])
      
    } else if(is.character(i)) {
      
      if(!all(i %in% names(x@mgnets))) {
        stop("No mgnet object with such a name in the list.")
      }
      return(x@mgnets[i])
      
    } else {
      stop("Index must be either numeric or a character vector of names.")
    }
  }
})