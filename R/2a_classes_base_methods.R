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
  if(length(object@abun!=0)) return(nrow(object@abun))
  else if(length(object@rela!=0)) return(nrow(object@rela))
  else if(length(object@norm!=0)) return(nrow(object@norm))
  else if(length(object@sample!=0)) return(nrow(object@sample))
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
  if(length(object@abun)!=0) return(ncol(object@abun))
  else if(length(object@taxa)!=0) return(nrow(object@taxa))
  else if(length(object@rela)!=0) return(ncol(object@rela))
  else if(length(object@norm)!=0) return(ncol(object@norm))
  else if(length(object@netw)!=0) return(vcount(object@netw))
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
  if(length(object@abun!=0)) return(rownames(object@abun))
  else if(length(object@sample!=0)) return(rownames(object@sample))
  else if(length(object@norm!=0)) return(rownames(object@norm))
  else if(length(object@rela!=0)) return(rownames(object@rela))
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
  if(length(object@abun!=0)) return(colnames(object@abun))
  else if(length(object@taxa)!=0) return(rownames(object@taxa))
  else if(length(object@rela)!=0) return(colnames(object@rela))
  else if(length(object@norm)!=0) return(colnames(object@norm))
  else if(length(object@netw)!=0) return(V(object@netw)$name)
  else return(character(length=0))
})

setMethod("taxa_id", "mgnetList", function(object) {
  sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE)
})


# INFO_SAMPLE_VARS
#------------------------------------------------------------------------------#
#' Get Sample Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `sample` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#' @export
#' @name sample_vars
#' @aliases sample_vars,mgnet-method sample_vars,mgnetList-method
setGeneric("sample_vars", function(object) standardGeneric("sample_vars"))

setMethod("sample_vars", "mgnet", function(object) {
  if(length(object@sample)!=0){
    return(colnames(object@sample))
  } else {
    return(character(length=0))
  }
})

setMethod("sample_vars", "mgnetList", function(object) {
  sapply(object@mgnets, sample_vars, simplify = FALSE, USE.NAMES = TRUE)
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
#' @aliases ncomm,mgnet-method ncomm,mgnetList-method
#' @export
setGeneric("ncomm", function(object) standardGeneric("ncomm"))

setMethod("ncomm", "mgnet", function(object){
  
  if(length(object@comm)!=0){
    
    sizes <- names(igraph::sizes(object@comm))
    class(sizes) <- "numeric"
    return(max(sizes))
    
  } else {
    
    return(numeric(0))
    
  }
  
})

setMethod("ncomm","mgnetList",
          function(object){
            sapply(object@mgnets, ncomm, simplify = TRUE, USE.NAMES = TRUE)
          })


# INFO_TAXA_VARS
#------------------------------------------------------------------------------#
#' Get Taxa Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `taxa` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#' @export
#' @name taxa_vars
#' @aliases taxa_vars,mgnet-method taxa_vars,mgnetList-method
setGeneric("taxa_vars", function(object) standardGeneric("taxa_vars"))

setMethod("taxa_vars", "mgnet", function(object) {
  if(length(object@taxa)!=0){
    return(colnames(object@taxa))
  } else {
    return(character(length=0))
  }
})

setMethod("taxa_vars", "mgnetList", function(object) {
  sapply(object@mgnets, taxa_vars, simplify = FALSE, USE.NAMES = TRUE)
})


# PULL INFO SAMPLE
#------------------------------------------------------------------------------#
#' Pull Sample Information from mgnet or mgnetList Objects
#'
#' @description
#' Retrieves specific sample information from the `sample` data frame within `mgnet` objects,
#' or each `mgnet` object within an `mgnetList`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param var The name or position of the column to retrieve from the `sample` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `mgnet` object, a vector containing the data from the specified column 
#'         of the `sample` data frame is returned. For an `mgnetList` object, a list of such vectors
#'         is returned, each corresponding to one `mgnet` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names.
#'
#' @export
#' @importFrom dplyr pull
#' @name pull_sample
#' @aliases pull_sample,mgnet-method pull_sample,mgnetList-method
setGeneric("pull_sample", function(object, var = -1) standardGeneric("pull_sample"))

setMethod("pull_sample", signature(object = "mgnet"), function(object, var = -1) {
  
  dplyr::pull(object@sample, {{var}})
})

setMethod("pull_sample", signature(object = "mgnetList"), function(object, var = -1) {
  sapply(object@mgnets, function(x) pull_sample(x, var), simplify = FALSE, USE.NAMES = TRUE)
})


# PULL INFO TAXA
#------------------------------------------------------------------------------#
#' Pull Taxa Information from mgnet or mgnetList Objects
#'
#' @description
#' Retrieves specific taxonomic information from the `taxa` data frame within `mgnet` objects,
#' or each `mgnet` object within an `mgnetList`. This function simplifies direct access to specific columns of interest
#' using dynamic column name handling.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param var The name or position of the column to retrieve from the `taxa` data frame. 
#'        If -1 (default), the last column of the data frame is returned. You can specify the column name unquoted due to non-standard evaluation.
#' @return For a single `mgnet` object, a vector containing the data from the specified column 
#'         of the `taxa` data frame is returned. For an `mgnetList` object, a list of such vectors
#'         is returned, each corresponding to one `mgnet` object in the list.
#'
#' @details
#' This function supports dynamic evaluation of the column name using `rlang` for unquoted names, allowing more
#' flexible and intuitive usage within data manipulation workflows.
#'
#' @export
#' @importFrom dplyr pull
#' @name pull_taxa
#' @aliases pull_taxa,mgnet-method pull_taxa,mgnetList-method
setGeneric("pull_taxa", function(object, var = -1) standardGeneric("pull_taxa"))

setMethod("pull_taxa", signature(object = "mgnet"), function(object, var = -1) {
  
  dplyr::pull(object@taxa, {{var}})
})

setMethod("pull_taxa", signature(object = "mgnetList"), function(object, var = -1) {
  sapply(object@mgnets, function(x) pull_taxa(x, var), simplify = FALSE, USE.NAMES = TRUE)
})


# comm_id
#------------------------------------------------------------------------------#
#' Retrieve Community IDs from mgnet or mgnetList Objects
#'
#' @description
#' This function retrieves the community IDs from an `mgnet` or `mgnetList` object. The community IDs are
#' derived from the `community` slot of the `mgnet` object or from each `mgnet` object within an `mgnetList`.
#' The function supports multiple output formats including list, data frame, and tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt The format of the output. Possible values are:
#' \itemize{
#'        \item `"list"`: Returns a named list where each element is a character vector of community IDs.
#'        \item `"df"`: Returns a data frame where each row corresponds to a taxa ID and each column to a community ID.
#'        \item `"tbl"`: Returns a tibble where each row corresponds to a taxa ID and each column to a community ID.
#'}
#'        Default is `"list"`.
#' @return For a single `mgnet` object, the function returns the community IDs in the specified format. 
#'         For an `mgnetList` object, it returns a list, data frame, or tibble, depending on the specified format, 
#'         where each entry corresponds to the community IDs for each `mgnet` object in the list.
#'         
#' @importFrom igraph membership
#' @importFrom stats setNames
#' @importFrom tibble rownames_to_column tibble
#' @aliases comm_id,mgnet-method comm_id,mgnetList-method
#' @export
setGeneric("comm_id", function(object, .fmt = "list") standardGeneric("comm_id"))

setMethod("comm_id", "mgnet", function(object, .fmt = "list"){
  
  if(length(object@comm) != 0){
    
    result <- switch(.fmt,
                     list = setNames(as.character(membership(object@comm)), taxa_id(object)),
                     df = data.frame("comm_id" = as.character(membership(object@comm)), row.names = taxa_id(object)),
                     tbl = tibble("taxa_id" = taxa_id(object), "comm_id" = as.character(membership(object@comm))))
    return(result)
    
  } else {
    
    return(character(0))
    
  }
})

setMethod("comm_id", "mgnetList", function(object, .fmt = "list"){
  .fmt <- match.arg(.fmt, choices = c("list", "df", "tbl"))
  
  if(.fmt == "list"){
    return(sapply(object@mgnets, function(x) comm_id(x), 
                  simplify = FALSE, USE.NAMES = TRUE))
  } else {
    taxa.merge <- unique(unlist(lapply(object@mgnets, taxa_id)))
    res <- matrix(NA_character_, nrow = length(taxa.merge), ncol = length(object@mgnets),
                  dimnames = list(taxa.merge, names(object@mgnets)))
    
    for(n in names(object@mgnets)){
      res[taxa_id(object@mgnets[[n]]), n] <- comm_id(object@mgnets[[n]])
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


# MGNETLIST ONLY METHODS
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Retrieve list of mgnet Objects from an mgnetList
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

#' Length of an mgnetList
#'
#' @description
#' Returns the number of `mgnet` objects contained within an `mgnetList` object.
#'
#' @param x An `mgnetList` object.
#' @return Integer value representing the number of `mgnet` objects in the `mgnetList`.
#' @export
setMethod("length", "mgnetList", function(x) {
  length(x@mgnets)
})


#' Names of mgnet Objects in an mgnetList
#'
#' @description
#' Retrieves the names of `mgnet` objects contained within an `mgnetList`.
#' Names provide a convenient way to reference and manage individual `mgnet` objects.
#'
#' @param x An `mgnetList` object.
#' @return A character vector of names of the `mgnet` objects.
#'         
#' @export
setMethod("names", "mgnetList", function(x) {
  names(x@mgnets)
})

#' Set Names of mgnet Objects in an mgnetList
#'
#' @description
#' Sets the names of `mgnet` objects contained within an `mgnetList`.
#' Names provide a convenient way to reference and manage individual `mgnet` objects.
#'
#' @param x An `mgnetList` object.
#' @param value A character vector representing the new names to be assigned to the `mgnet` objects.
#' @return The modified `mgnetList` object with updated names.
#'         
#' @importFrom methods validObject
#' @export
setMethod("names<-", "mgnetList", function(x, value) {
  names(x@mgnets) <- value
  validObject(x)
  x
})