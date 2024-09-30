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
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods new
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix
#' @export
setMethod(f="[", signature="mgnet",function(x, i, j) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Helper function to check uniqueness (only for numeric or character indices)
  check_unique <- function(indices, label) {
    if (any(duplicated(indices))) {
      stop(sprintf("Error: '%s' contains duplicate values.", label))
    }
  }
  
  # Validate and process 'i' indices (samples)
  if (!missing(i)) {
    if (is.numeric(i)) {
      if (any(i < 1 | i > nsample(x))) {
        stop("Error: Sample indices 'i' out of range.")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.character(i)) {
      i <- match(i, sample_id(x))
      if (any(is.na(i))) {
        stop("Error: Some sample names in 'i' not found in sample_id.")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.logical(i)) {
      if (length(i) != nsample(x)) {
        stop("Error: Logical index 'i' must have the same length as number of samples.")
      }
      i <- which(i)  # Convert logical to numeric
    } else {
      stop("Error: Invalid sample indices 'i'. Must be numeric, character, or logical.")
    }
  } else {
    i <- seq_len(nsample(x))  # Default: all samples
  }
  
  # Validate and process 'j' indices (taxa)
  if (!missing(j)) {
    if (is.numeric(j)) {
      if (any(j < 1 | j > ntaxa(x))) {
        stop("Error: Taxa indices 'j' out of range.")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.character(j)) {
      j <- match(j, taxa_id(x))
      if (any(is.na(j))) {
        stop("Error: Some taxa names in 'j' not found in taxa_id.")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.logical(j)) {
      if (length(j) != ntaxa(x)) {
        stop("Error: Logical index 'j' must have the same length as number of taxa.")
      }
      j <- which(j)  # Convert logical to numeric
    } else {
      stop("Error: Invalid taxa indices 'j'. Must be numeric, character, or logical.")
    }
  } else {
    j <- seq_len(ntaxa(x))  # Default: all taxa
  }
  
  # Determine if full sample set is used to retain network and community
  full_i <- all(seq_len(nsample(x)) %in% i) || all(sample_id(x) %in% i)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # Handle empty indices cases
  if ( length(i) == 0 && length(j) == 0 ){
    
    return(new("mgnet")) # Return a completely new, empty mgnet object
  
  } else if ( length(i) == 0 && length(j) > 0 ){
    
    # Case where only taxa are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(nrow=0, ncol=0)
    new_mgnet@rela <- matrix(nrow=0, ncol=0)
    new_mgnet@norm <- matrix(nrow=0, ncol=0)
    new_mgnet@meta <- data.frame()
    new_mgnet@taxa <- x@taxa[j,]
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    validObject(new_mgnet)
    return(new_mgnet)
    
  } else if ( length(i) > 0 && length(j) == 0 ){
    
    # Case where only samples are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@rela <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@norm <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@meta <- x@meta[i,]
    new_mgnet@taxa <- data.frame()
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    validObject(new_mgnet)
    return(new_mgnet)
    
  } else {
    
    # Standard case where both i and j are present
    # Initialize new data subsets
    if(length(x@abun) != 0) abun.new <- x@abun[i,j,drop=F] else abun.new <- x@abun
    if(length(x@rela) != 0) rela.new <- x@rela[i,j,drop=F] else rela.new <- x@rela
    if(length(x@norm) != 0) norm.new <- x@norm[i,j,drop=F] else norm.new <- x@norm
    if(length(x@meta) != 0) meta.new <- x@meta[i, ,drop=F] else meta.new <- x@meta
    if(length(x@taxa) != 0) taxa.new <- x@taxa[j, ,drop=F] else taxa.new <- x@taxa

    if(length(x@netw) == 0){
      
      # When network is missing only abundances and metadata are returned
      return(mgnet(abun = abun.new,
                   rela = rela.new,
                   norm = norm.new,
                   meta = meta.new,
                   taxa = taxa.new))
      
    } else if (length(x@netw)!=0 & full_i) {
      
      # Obtain the adjacency matrix
      if(igraph::is_weighted(x@netw)){
        adj <- igraph::as_adjacency_matrix(x@netw, attr = "weight", sparse = FALSE)
      } else {
        adj <- igraph::as_adjacency_matrix(x@netw, sparse = FALSE)
      }
      
      #Reorder vertices
      adj_reorder <- adj[j,j]
      
      # Reconstruct the graph from the reordered adjacency matrix
      netw.new <- igraph::graph_from_adjacency_matrix(adj_reorder, mode = "undirected", weighted = TRUE)
      
      
      if(length(x@comm) != 0){
        comm.new <- x@comm
        if(is.character(j)) j <- which(taxa_id(x)%in%j)
        comm.new$membership <- x@comm$membership[j]
        comm.new$vcount <- length(comm.new$membership)
        comm.new$modularity <- NA
      } else {
        comm.new <- x@comm
      }
      
      return(mgnet(abun = abun.new,
                   rela = rela.new,
                   norm = norm.new,
                   meta = meta.new,
                   taxa = taxa.new,
                   netw = netw.new,
                   comm = comm.new))
      
    } else if (length(x@netw)!=0 & !full_i) {
      # SUBCASE NETWORK PRESENT AND I PRESENT
      
      warning("Subsetting by samples removes 'netw' and 'comm' slots.\n")
      return(mgnet(abun = abun.new,
                   rela = rela.new,
                   norm = norm.new,
                   meta = meta.new,
                   taxa = taxa.new))
    
    } else {
      
      stop("why are you here? only for suffering? (cit. Kaz).")
      
    }
    
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


# MGNETS
#------------------------------------------------------------------------------#
#' Retrieve mgnet Objects from an mgnetList
#'
#' This method extracts and returns the list of `mgnet` objects contained within an `mgnetList` object,
#' offering direct access to the individual `mgnet` objects for further analysis or manipulation. It functions
#' similarly to the \code{\link{as.list}} method for `mgnetList` objects, both providing a way to access 
#' the contained `mgnet` objects.
#'
#' @param object An `mgnetList` object from which `mgnet` objects are to be retrieved.
#' 
#' @return A list of `mgnet` objects, where each element is an individual `mgnet` object. The elements
#'         are named, corresponding to the names of the `mgnet` objects within the `mgnetList`, allowing
#'         for easy identification and access. This facilitates individual or batch processing of 
#'         metagenomic network data.
#'
#' @seealso
#' \code{\link[=mgnet]{mgnet}} for details on the `mgnet` class and functionalities.
#' \code{\link[=mgnetList]{mgnetList}} for guidance on creating and managing `mgnetList` objects.
#' \code{\link{as.list}} for a similar method that converts an `mgnetList` into a list of `mgnet` objects,
#' serving the same purpose as `mgnets`.
#'
#' @export
#' @name mgnets
#' @aliases mgnets,mgnetList-method
setGeneric("mgnets", function(object) standardGeneric("mgnets"))

setMethod("mgnets", "mgnetList", function(object) {
  object@mgnets
})