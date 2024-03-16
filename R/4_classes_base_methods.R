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
  else if(length(object@abundance!=0)) return(nrow(object@abundance))
  else if(length(object@info_sample!=0)) return(nrow(object@info_sample))
  else return(NULL)
})

setMethod("nsample", "mgnetList", function(object) {
  sapply(object@mgnetObjects, nsample, simplify=T, USE.NAMES=T)
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
  else if(length(object@log_abundance)!=0) return(ncol(object@log_abundance))
  else if(length(object@network)!=0) return(vcount(object@network))
  else return(NULL)
})

setMethod("ntaxa", "mgnetList", function(object) {
  sapply(object@mgnetObjects, ntaxa, simplify=T, USE.NAMES=T)
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
  else if(length(object@log_abundance!=0)) return(rownames(object@log_abundance))
  else return(NULL)
})

setMethod("sample_id", "mgnetList", function(object) {
  sapply(object@mgnetObjects, sample_id, simplify=F, USE.NAMES=T)
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
  if(length(object@data!=0)) return(colnames(object@data))
  else if(length(object@taxa!=0)) return(rownames(object@taxa))
  else if(length(object@meta_taxa)!=0) return(rownames(object@meta_taxa))
  else if(length(object@log_data)!=0) return(colnames(object@log_data))
  else if(length(object@netw)!=0) return(V(object@netw)$name)
  else return(NULL)
})

setMethod("taxa_id", "mgnetList", function(object) {
  sapply(object@mgnetObjects, taxa_id, simplify=F, USE.NAMES=T)
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
  colnames(object@lineage)
})

setMethod("ranks", "mgnetList", function(object) {
  sapply(object@mgnetObjects, ranks, simplify=T, USE.NAMES=T)
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
  colnames(object@info_sample)
})

setMethod("info_sample_vars", "mgnetList", function(object) {
  sapply(object@mgnetObjects, info_sample_vars, simplify=F, USE.NAMES=T)
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
  colnames(object@info_taxa)
})

setMethod("info_taxa_vars", "mgnetList", function(object) {
  sapply(object@mgnetObjects, info_taxa_vars, simplify=F, USE.NAMES=T)
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
#' @aliases taxa_name,mgnet-method taxa_name,mgnetList-method
setGeneric("taxa_name", function(object, rank = "missing") standardGeneric("taxa_name"))

#' @rdname taxa_name
setMethod("taxa_name", signature(object = "mgnet", rank = "missing"), function(object) {
  finest_rank <- colnames(object@lineage)[ncol(object@lineage)]
  object@lineage[, finest_rank]
})

#' @rdname taxa_name
setMethod("taxa_name", signature(object = "mgnet", rank = "character"), function(object, rank) {
  if(!(rank %in% colnames(object@lineage))) {
    stop("Specified rank is not available. Available ranks are: ", toString(colnames(object@lineage)))
  }
  object@lineage[, rank]
})

#' @rdname taxa_name
setMethod("taxa_name", signature(object = "mgnetList", rank = "missing"), function(object) {
  sapply(object@mgnetObjects, function(x) taxa_name(x), simplify=F, USE.NAMES=T)
})

#' @rdname taxa_name
setMethod("taxa_name", signature(object = "mgnetList", rank = "character"), function(object, rank) {
  sapply(object@mgnetObjects, function(x) taxa_name(x, rank), simplify=F, USE.NAMES=T)
})
