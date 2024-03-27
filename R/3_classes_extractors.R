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

#' Access mgnet Objects by Name
#'
#' This method allows for easy access to individual `mgnet` objects within an `mgnetList` by their names.
#'
#' @param x An `mgnetList` object.
#' @param name The name of the `mgnet` object to access.
#' @return The `mgnet` object with the specified name.
#' 
#' @importFrom methods slot
#' @export
setMethod("$", "mgnetList", function(x, name) {
  slot(x, "mgnets")[[name]]
})


#' Subset mgnetList
#'
#' This method subsets an `mgnetList` object, returning a new `mgnetList` containing only the selected `mgnet` objects.
#'
#' @param x An `mgnetList` object.
#' @param i Indices or names specifying the `mgnet` objects to retain.
#' @param j Not used.
#' @param ... Additional arguments (currently unused).
#' @param drop Not used.
#' @return A new `mgnetList` object containing the selected `mgnet` objects.
#' 
#' @importFrom methods slot
#' @export
setMethod("[", "mgnetList", function(x, i, j, ..., drop = FALSE) {
  
  if (!missing(...)) {
    stop("Only the 'i' argument is supported for extracting from mgnetList objects with `[`.")
  }
  
  if (!missing(j)) {
    stop("The 'j' argument is not supported for mgnetList objects.")
  }
  if (!isFALSE(drop)) {
    stop("The 'drop' argument is not supported for mgnetList objects. mgnetList does not support dropping dimensions.")
  }
  
  new_mgnets <- slot(x, "mgnets")[i]
  new("mgnetList", mgnets=new_mgnets)
})


#' Extract a Single mgnet Object
#'
#' This method extracts a single `mgnet` object from an `mgnetList` by index or name.
#'
#' @param x An `mgnetList` object.
#' @param i Index or name specifying the `mgnet` object to extract.
#' @param ... Additional arguments (currently unused).
#' @param exact Logical indicating if exact matching is required (default is `TRUE`).
#' @return The extracted `mgnet` object.
#' 
#' @importFrom methods slot
#' @export
setMethod("[[", "mgnetList", function(x, i, ..., exact = TRUE) {
  
  if (!missing(...) || !isTRUE(exact)) {
    stop("Only the 'i' argument is supported for extracting from mgnetList objects with `[[`. 'exact' must be TRUE.")
  }
  
  slot(x, "mgnets")[[i, exact = exact]]
})


#' Replace a Single mgnet Object
#'
#' This method replaces a single `mgnet` object within an `mgnetList` by index or name.
#'
#' @param x An `mgnetList` object.
#' @param i Index or name specifying which `mgnet` object to replace.
#' @param value The new `mgnet` object to be inserted into the list.
#' @return The `mgnetList` object with the specified `mgnet` object replaced.
#' 
#' @importFrom methods slot<-
#' @export
setMethod("[[<-", "mgnetList", function(x, i, value) {
  if (!inherits(value, "mgnet")) {
    stop("The 'value' must be an instance of the 'mgnet' class.")
  }
  
  # Extract the current list of mgnet objects
  mgnets <- slot(x, "mgnets")
  
  # Replace the specified mgnet object
  mgnets[[i]] <- value
  
  # Update the slot with the modified list
  slot(x, "mgnets") <- mgnets
  
  # Validate the modified mgnetList
  validObject(x)
  
  return(x)
})
