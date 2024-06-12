# NSAMPLE
#------------------------------------------------------------------------------#
#' Get Number of Samples
#'
#' Returns an integer indicating the number of samples in an `mgnet` object
#' or a list of integers for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, an integer representing the number of samples.
#'         For an `mgnetList` object, a list of integers, each representing the number
#'         of samples in the corresponding `mgnet` objects.
#' @export
#' @name nsample
#' @aliases nsample,mgnet-method nsample,mgnetList-method
setGeneric("nsample", function(object) standardGeneric("nsample"))

setMethod("nsample", "mgnet", function(object) {
  if(length(object@abundance!=0)) return(nrow(object@abundance))
  else if(length(object@rel_abundance!=0)) return(nrow(object@rel_abundance))
  else if(length(object@norm_abundance!=0)) return(nrow(object@norm_abundance))
  else if(length(object@info_sample!=0)) return(nrow(object@info_sample))
  else return(0)
})

setMethod("nsample", "mgnetList", function(object) {
  sapply(object@mgnets, nsample, simplify = TRUE, USE.NAMES = TRUE)
})


# NTAXA
#------------------------------------------------------------------------------#
#' Get Number of Taxa
#'
#' Returns an integer indicating the number of taxa in an `mgnet` object
#' or a list of integers for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, an integer representing the number of taxa.
#'         For an `mgnetList` object, a list of integers, each representing the number
#'         of taxa in the corresponding `mgnet` objects.
#' @export
#' @importFrom igraph vcount
#' @name ntaxa
#' @aliases ntaxa,mgnet-method ntaxa,mgnetList-method
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))

setMethod("ntaxa", "mgnet", function(object) {
  if(length(object@abundance)!=0) return(ncol(object@abundance))
  else if(length(object@lineage)!=0) return(nrow(object@lineage))
  else if(length(object@info_taxa)!=0) return(nrow(object@info_taxa))
  else if(length(object@rel_abundance)!=0) return(ncol(object@rel_abundance))
  else if(length(object@norm_abundance)!=0) return(ncol(object@norm_abundance))
  else if(length(object@network)!=0) return(vcount(object@network))
  else return(0)
})

setMethod("ntaxa", "mgnetList", function(object) {
  sapply(object@mgnets, ntaxa, simplify = TRUE, USE.NAMES = TRUE)
})


# SAMPLE_ID
#------------------------------------------------------------------------------#
#' Get Sample IDs
#'
#' Returns the names of samples as a character vector from an `mgnet` object
#' or lists of sample IDs for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector representing the IDs of samples.
#'         For an `mgnetList` object, a list of character vectors, each representing the sample IDs
#'         in the corresponding `mgnet` objects.
#' @export
#' @name sample_id
#' @aliases sample_id,mgnet-method sample_id,mgnetList-method
setGeneric("sample_id", function(object) standardGeneric("sample_id"))

setMethod("sample_id", "mgnet", function(object) {
  if(length(object@abundance!=0)) return(rownames(object@abundance))
  else if(length(object@info_sample!=0)) return(rownames(object@info_sample))
  else if(length(object@norm_abundance!=0)) return(rownames(object@norm_abundance))
  else if(length(object@rel_abundance!=0)) return(rownames(object@rel_abundance))
  else return(character(length=0))
})

setMethod("sample_id", "mgnetList", function(object) {
  sapply(object@mgnets, sample_id, simplify = FALSE, USE.NAMES = TRUE)
})


# TAXA_ID
#------------------------------------------------------------------------------#
#' Get Taxa IDs
#'
#' Retrieves the IDs of taxa from an `mgnet` object or lists of taxa IDs for each 
#' `mgnet` object within an `mgnetList`. Taxa IDs represent unique identifiers 
#' for the taxa present in the dataset.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector representing the IDs of taxa.
#'         For an `mgnetList` object, a list of character vectors, each representing 
#'         the taxa IDs in the corresponding `mgnet` objects.
#' @export
#' @importFrom igraph V
#' @name taxa_id
#' @aliases taxa_id,mgnet-method taxa_id,mgnetList-method
setGeneric("taxa_id", function(object) standardGeneric("taxa_id"))

setMethod("taxa_id", "mgnet", function(object) {
  if(length(object@abundance!=0)) return(colnames(object@abundance))
  else if(length(object@lineage!=0)) return(rownames(object@lineage))
  else if(length(object@info_taxa)!=0) return(rownames(object@info_taxa))
  else if(length(object@rel_abundance)!=0) return(colnames(object@rel_abundance))
  else if(length(object@norm_abundance)!=0) return(colnames(object@norm_abundance))
  else if(length(object@network)!=0) return(V(object@network)$name)
  else return(character(length=0))
})

setMethod("taxa_id", "mgnetList", function(object) {
  sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE)
})


# RANKS
#------------------------------------------------------------------------------#
#' Get Taxonomic Ranks
#'
#' Retrieves the taxonomic classification ranks from an `mgnet` object or lists 
#' of taxonomic ranks for each `mgnet` object within an `mgnetList`. 
#' Taxonomic ranks are derived from the column names of the `lineage` slot, 
#' representing different levels of taxonomic classification.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of taxonomic ranks.
#'         For an `mgnetList` object, a list of character vectors, each representing 
#'         the taxonomic ranks in the corresponding `mgnet` objects.
#' @export
#' @name ranks
#' @aliases ranks,mgnet-method ranks,mgnetList-method
setGeneric("ranks", function(object) standardGeneric("ranks"))

setMethod("ranks", "mgnet", function(object) {
  
  if(length(object@lineage)!=0) return(colnames(object@lineage))
  else return(character(length=0))
  
})

setMethod("ranks", "mgnetList", function(object) {
  sapply(object@mgnets, ranks, simplify = FALSE, USE.NAMES = TRUE)
})


# INFO_SAMPLE_VARS
#------------------------------------------------------------------------------#
#' Get Sample Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `info_sample` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#' @export
#' @name info_sample_vars
#' @aliases info_sample_vars,mgnet-method info_sample_vars,mgnetList-method
setGeneric("info_sample_vars", function(object) standardGeneric("info_sample_vars"))

setMethod("info_sample_vars", "mgnet", function(object) {
  if(length(info_sample)!=0){
    return(colnames(object@info_sample))
  } else {
    return(character(length=0))
  }
})

setMethod("info_sample_vars", "mgnetList", function(object) {
  sapply(object@mgnets, info_sample_vars, simplify = FALSE, USE.NAMES = TRUE)
})


# INFO_TAXA_VARS
#------------------------------------------------------------------------------#
#' Get Taxa Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `info_taxa` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#' @export
#' @name info_taxa_vars
#' @aliases info_taxa_vars,mgnet-method info_taxa_vars,mgnetList-method
setGeneric("info_taxa_vars", function(object) standardGeneric("info_taxa_vars"))

setMethod("info_taxa_vars", "mgnet", function(object) {
  if(length(info_taxa)!=0){
    return(colnames(object@info_taxa))
  } else {
    return(character(length=0))
  }
})

setMethod("info_taxa_vars", "mgnetList", function(object) {
  sapply(object@mgnets, info_taxa_vars, simplify = FALSE, USE.NAMES = TRUE)
})


# TAXA_NAME
#------------------------------------------------------------------------------#
#' Retrieve Taxa Names at Specified Rank
#'
#' Extracts the names of taxa at the specified taxonomic rank from an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`. Unlike `taxa_id` which returns
#' identifiers associated with the row/column names of the input matrix/data.frame
#' (often reflecting the finest classification, such as OTUs, where IDs and names may coincide),
#' `taxa_name` delves into the `lineage` matrix to fetch descriptive names for the entities at the
#' chosen rank. This provides more characteristic names of the entities.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param rank The taxonomic rank for which names are to be retrieved. If not specified,
#'        the function returns names at the finest available taxonomic rank. For `mgnet` objects,
#'        this is derived from the column names of the `lineage` matrix. For `mgnetList` objects,
#'        it iterates over each contained `mgnet` object to perform the extraction.
#' @return For an `mgnet` object, a character vector of taxa names at the specified rank.
#'         For an `mgnetList` object, a list where each element is a character vector
#'         of taxa names at the specified rank from each `mgnet` object within the list.
#' @export
#' @exportMethod taxa_name
#' @name taxa_name
#' @aliases taxa_name,mgnet,missing-method taxa_name,mgnet,character-method taxa_name,mgnetList,missing-method taxa_name,mgnetList,character-method
setGeneric("taxa_name", function(object, rank = "missing") standardGeneric("taxa_name"))

setMethod("taxa_name", signature(object = "mgnet", rank = "missing"), function(object) {
  
  if(length(object@lineage) == 0 ) return(character(length = 0))
  return(object@lineage[, ncol(lineage)])
  
})

setMethod("taxa_name", signature(object = "mgnet", rank = "character"), function(object, rank) {
  
  if(length(object@lineage)==0) return(character(length = 0))
  
  if(!(rank %in% colnames(object@lineage))) {
    stop("Specified rank is not available. Available ranks are: ", toString(colnames(object@lineage)))
  }
  object@lineage[, rank]
})

setMethod("taxa_name", signature(object = "mgnetList", rank = "missing"), function(object) {
  sapply(object@mgnets, function(x) taxa_name(x), simplify=F, USE.NAMES=T)
})

setMethod("taxa_name", signature(object = "mgnetList", rank = "character"), function(object, rank) {
  sapply(object@mgnets, function(x) taxa_name(x, rank), simplify=F, USE.NAMES=T)
})


# community_id
#------------------------------------------------------------------------------#
#' Retrieve Community Memberships for Network Vertices
#'
#' This function fetches the community membership IDs for each vertex in a network or a list of networks. 
#' It supports `mgnet` objects, representing single networks, and `mgnetList` objects, representing collections of networks. 
#' The function can return the results in various formats, including lists, data frames, or tibbles, depending on user preference.
#'
#' @param object An object of class `mgnet` or `mgnetList`. For `mgnet`, the function returns the community membership IDs 
#' of all vertices in the network. For `mgnetList`, the function operates on each network in the list and compiles the results.
#' @param .fmt A character string indicating the desired output format for the community membership IDs when the input 
#' object is an `mgnetList`. Options include "list" for a list of vectors, "df" for a data frame, and "tbl" for a tibble. 
#' The default is "list". This parameter is ignored if the input is an `mgnet` object.
#'
#' @return Depending on the input object and the `.fmt` parameter:
#' \itemize{
#'   \item For an `mgnet` object, returns a named vector where each name is a vertex ID and each value is the corresponding community membership ID.
#'   \item For an `mgnetList` object with `.fmt` set to "list", returns a list of named vectors, each corresponding to one network in the list.
#'   \item For an `mgnetList` object with `.fmt` set to "df", returns a data frame where each row corresponds to a vertex (across all networks) and each column corresponds to one network.
#'   \item For an `mgnetList` object with `.fmt` set to "tbl", returns a tibble similar to the data frame format, but with `taxa_id` as the first column.
#' }
#'
#' @importFrom igraph membership
#' @importFrom stats setNames
#' @importFrom tibble rownames_to_column tibble
#' @aliases community_members,mgnet-method community_members,mgnetList-method
#' @export
#'
#' @details The function leverages the `membership` function from the `igraph` package to determine the community membership of vertices.
#' For `mgnetList` objects, the function applies the operation to each network in the list and organizes the results according to the specified format.
#' When dealing with `mgnetList` objects and choosing "df" or "tbl" for `.fmt`, the function uniquely identifies each vertex across all networks to compile the community memberships coherently.
#' 
setGeneric("community_members", function(object, .fmt = "list") standardGeneric("community_members"))

setMethod("community_members", "mgnet", function(object, .fmt = "list"){
  
  if(length(object@community) != 0){
    
    result <- switch(.fmt,
                     list = setNames(as.character(membership(object@community)), taxa_id(object)),
                     df = data.frame("community" = as.character(membership(object@community)), row.names = taxa_id(object)),
                     tbl = tibble("taxa_id" = taxa_id(object), "community" = as.character(membership(object@community))))
    return(result)
    
  } else {
    
    return(character(0))
    
  }
  
})

setMethod("community_members", "mgnetList", function(object, .fmt = "list"){
  .fmt <- match.arg(.fmt, choices = c("list", "df", "tbl"))
  
  if(.fmt == "list"){
    return(sapply(object@mgnets, function(x) community_members(x), 
                  simplify = FALSE, USE.NAMES = TRUE))
  } else {
    taxa.merge <- unique(unlist(lapply(object@mgnets, taxa_id)))
    res <- matrix(NA_character_, nrow = length(taxa.merge), ncol = length(object@mgnets),
                  dimnames = list(taxa.merge, names(object@mgnets)))
    
    for(n in names(object@mgnets)){
      res[taxa_id(object@mgnets[[n]]), n] <- community_members(object@mgnets[[n]])
    }
    
    if(.fmt == "df"){
      return(data.frame(res, stringsAsFactors = FALSE))
    } else if(.fmt == "tbl"){
      # Convert to tibble with taxa_id as the first column
      res_df <- data.frame(res, stringsAsFactors = FALSE)
      res_tbl <- tibble::rownames_to_column(res_df, var = "taxa_id")
      return(res_tbl)
    }
  }
})


# NCOMMUNITY
#------------------------------------------------------------------------------#
#' Get the Number of Communities
#'
#' Retrieves the number of communities present in the `comm` slot of an `mgnet` object or each `mgnet` object within an `mgnetList`.
#'
#' @description
#' This function returns an integer indicating the number of communities detected in a network analysis represented by an `mgnet` object, 
#' or a vector of integers for an `mgnetList` object, where each element corresponds to the number of communities in each respective `mgnet` object.
#' For `mgnet` objects without any communities, the function returns 0. For `mgnetList` objects, it iteratively applies this logic to each contained `mgnet` object.
#'
#' @param object An object of class `mgnet` or `mgnetList`. For `mgnet`, it calculates the number of communities directly from the `comm` slot.
#' For `mgnetList`, it applies the calculation to each `mgnet` object within the list and returns a vector of results.
#'
#' @return For an `mgnet` object, returns an integer representing the number of communities.
#' For an `mgnetList` object, returns a numeric vector with each element representing the number of communities in each `mgnet` object.
#'
#' @export
#' @importFrom igraph sizes
#' @aliases ncommunity,mgnet-method ncommunity,mgnetList-method
#' @export
setGeneric("ncommunity", function(object) standardGeneric("ncommunity"))

setMethod("ncommunity", "mgnet", function(object){
  
  if(length(object@community)!=0){
    
    sizes <- names(igraph::sizes(object@community))
    class(sizes) <- "numeric"
    return(max(sizes))
    
  } else {
    
    return(numeric(0))
    
  }
  
})

setMethod("ncommunity","mgnetList",
          function(object){
            sapply(object@mgnets, ncommunity, simplify = TRUE, USE.NAMES = TRUE)
          })


# SAMPLE SUM
#------------------------------------------------------------------------------#
#' Calculate Sample Sum of Abundance Data in mgnet Objects
#'
#' This function calculates the sum of abundance values for each sample within the
#' abundance matrix of an `mgnet` object. It returns a numeric vector where each element
#' corresponds to the sum of abundance values for a sample.
#'
#' @param object An `mgnet` object containing an abundance matrix.
#' @param na.rm Logical indicating whether NA values should be removed before 
#'        summing the abundance values. Defaults to FALSE.
#' @return A numeric vector with the sum of abundance values for each sample.
#' @export
#' @aliases sample_sum,mgnet-method sample_sum,mgnetList-method
#' @export
setGeneric("sample_sum", function(object, na.rm = FALSE) standardGeneric("sample_sum"))
setMethod("sample_sum", "mgnet", function(object, na.rm = FALSE) {
  if(is.null(object@abundance) || nrow(object@abundance) == 0) {
    stop("The abundance matrix is missing or empty.")
  }
  
  # Calculate the sum of abundance values for each sample
  sample_sums <- rowSums(object@abundance, na.rm = na.rm)
  
  return(sample_sums)
})

setMethod("sample_sum","mgnetList",
          function(object, na.rm = FALSE){
            sapply(object@mgnets, sample_sum, simplify = TRUE, USE.NAMES = TRUE)
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