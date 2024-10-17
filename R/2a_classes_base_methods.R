# NSAMPLE
#------------------------------------------------------------------------------#
#' Get Number of Samples
#'
#' Returns an integer indicating the number of samples in an `mgnet` object
#' or a named vector of integers for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, an integer representing the number of samples.
#'         For an `mgnetList` object, a named vector of integers, each representing the number
#'         of samples in the corresponding `mgnet` objects.
#'         
#' @examples
#' data(HMP2, package = "mgnet")
#' nsample(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' nsample(subjects_HMP2)  
#' 
#' @export
#' @name nsample
#' @aliases nsample,mgnet-method nsample,mgnetList-method
setGeneric("nsample", function(object) standardGeneric("nsample"))

setMethod("nsample", "mgnet", function(object) {
  if(length(object@abun!=0)) return(nrow(object@abun))
  else if(length(object@rela!=0)) return(nrow(object@rela))
  else if(length(object@norm!=0)) return(nrow(object@norm))
  else if(length(object@meta!=0)) return(nrow(object@meta))
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
#' or a named vector of integers for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, an integer representing the number of taxa.
#'         For an `mgnetList` object, a named vector of integers, each representing the number
#'         of taxa in the corresponding `mgnet` objects.
#'         
#' @examples
#' data(HMP2, package = "mgnet")
#' ntaxa(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' ntaxa(subjects_HMP2)  
#' 
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
#' Get sample IDs from mgnet or mgnetList Objects
#'
#' Retrieves the IDs of sample from an `mgnet` object or lists of sample IDs for each 
#' `mgnet` object within an `mgnetList`. For `mgnetList` objects, this method can 
#' output the results either as a list or as a organized tibble.
#'
#' @param object An `mgnet` or `mgnetList` object from which to extract sample IDs.
#' @param .fmt Character; specifies the output format when retrieving sample IDs from 
#'        an `mgnetList`. Accepted values are \code{"list"} for a list of character vectors, 
#'        and \code{"tbl"} for a tibble. The default is \code{"list"}. The \code{"tbl"} format
#'        returns a tibble with two columns: \code{"mgnet"}, indicating the name of the `mgnet` object
#'        each sample ID originates from, and \code{"sample_id"}, listing the sample IDs themselves.
#'        
#' @importFrom purrr imap
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
#' @name sample_id
#' @aliases sample_id,mgnet-method sample_id,mgnetList-method
#' @examples
#' data(HMP2, package = "mgnet")
#' sample_id(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' sample_id(subjects_HMP2, .fmt = "list")  
#' sample_id(subjects_HMP2, .fmt = "tbl")  
#' 
#' @seealso \link{mgnet} and \link{mgnetList} for details on the classes.
setGeneric("sample_id", function(object, .fmt = "list") standardGeneric("sample_id"))

setMethod("sample_id", "mgnet", function(object, .fmt = "list") {
  if(length(object@abun) != 0) return(rownames(object@abun))
  else if(length(object@meta) != 0) return(rownames(object@meta))
  else if(length(object@norm) != 0) return(rownames(object@norm))
  else if(length(object@rela) != 0) return(rownames(object@rela))
  else return(character(0))
})

setMethod("sample_id", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "tbl"))
  
  if (.fmt == "list") {
    
    return(sapply(object@mgnets, sample_id, simplify = FALSE, USE.NAMES = TRUE))
    
  } else {
    
    sapply(object@mgnets, sample_id, simplify = FALSE, USE.NAMES = TRUE) %>%
      purrr::imap(\(x,y) tibble::tibble("mgnet" = y, "sample_id" = x)) %>%
      dplyr::bind_rows() %>%
      return()
    
  }
})



# TAXA_ID
#------------------------------------------------------------------------------#
#' Get Taxa IDs from mgnet or mgnetList Objects
#'
#' Retrieves the IDs of taxa from an `mgnet` object or lists of taxa IDs for each 
#' `mgnet` object within an `mgnetList`. For `mgnetList` objects, this method can 
#' output the results either as a list or as a organized tibble.
#'
#' @param object An `mgnet` or `mgnetList` object from which to extract taxa IDs.
#' @param .fmt Character; specifies the output format when retrieving taxa IDs from 
#'        an `mgnetList`. Accepted values are \code{"list"} for a list of character vectors, 
#'        and \code{"tbl"} for a tibble. The default is \code{"list"}. The \code{"tbl"} format
#'        returns a tibble with two columns: \code{"mgnet"}, indicating the name of the `mgnet` object
#'        each taxa ID originates from, and \code{"taxa_id"}, listing the taxa IDs themselves.
#'        
#' @importFrom purrr imap
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @export
#' @name taxa_id
#' @aliases taxa_id,mgnet-method taxa_id,mgnetList-method
#' @examples
#' data(HMP2, package = "mgnet")
#' taxa_id(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' taxa_id(subjects_HMP2, .fmt = "list")  
#' taxa_id(subjects_HMP2, .fmt = "tbl")  
#' 
#' @seealso \link{mgnet} and \link{mgnetList} for details on the classes.
setGeneric("taxa_id", function(object, .fmt = "list") standardGeneric("taxa_id"))

setMethod("taxa_id", "mgnet", function(object, .fmt = "list") {
  if(length(object@abun) != 0) return(colnames(object@abun))
  else if(length(object@taxa) != 0) return(rownames(object@taxa))
  else if(length(object@norm) != 0) return(colnames(object@norm))
  else if(length(object@rela) != 0) return(colnames(object@rela))
  else return(character(0))
})

setMethod("taxa_id", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "tbl"))
  
  if (.fmt == "list") {
    
    return(sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE))
    
  } else {
    
    sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE) %>%
      purrr::imap(\(x,y) tibble::tibble("mgnet" = y, "taxa_id" = x)) %>%
      dplyr::bind_rows() %>%
      return()
    
  }
})

# INFO_SAMPLE_VARS
#------------------------------------------------------------------------------#
#' Get Sample Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `meta` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#'         
#' @examples
#' data(HMP2, package = "mgnet")
#' meta_vars(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' meta_vars(subjects_HMP2)  
#' 
#' @export
#' @name meta_vars
#' @aliases meta_vars,mgnet-method meta_vars,mgnetList-method
setGeneric("meta_vars", function(object) standardGeneric("meta_vars"))

setMethod("meta_vars", "mgnet", function(object) {
  if(length(object@meta)!=0){
    return(colnames(object@meta))
  } else {
    return(character(length=0))
  }
})

setMethod("meta_vars", "mgnetList", function(object) {
  sapply(object@mgnets, meta_vars, simplify = FALSE, USE.NAMES = TRUE)
})


# NCOMMUNITY
#------------------------------------------------------------------------------#
#' Get the Number of Communities
#'
#' Retrieves the number of communities present in the `comm` slot of an `mgnet` object or each `mgnet` object within an `mgnetList`.
#'
#' @description
#' This function returns an integer indicating the number of communities detected in a network analysis represented by an `mgnet` object, 
#' or a named vector of integers for an `mgnetList` object, where each element corresponds to the number of communities in each respective `mgnet` object.
#' For `mgnet` objects without any communities, the function returns 0. For `mgnetList` objects, it iteratively applies this logic to each contained `mgnet` object.
#'
#' @param object An object of class `mgnet` or `mgnetList`. For `mgnet`, it calculates the number of communities directly from the `comm` slot.
#' For `mgnetList`, it applies the calculation to each `mgnet` object within the list and returns a vector of results.
#'
#' @return For an `mgnet` object, returns an integer representing the number of communities.
#' For an `mgnetList` object, returns a numeric named vector with each element representing the number of communities in each `mgnet` object.
#' 
#' @examples
#' data(HMP2, package = "mgnet")
#' ncomm(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' ncomm(subjects_HMP2)  
#' 
#' @export
#' @name ncomm
#' @aliases ncomm,mgnet-method ncomm,mgnetList-method
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
#'         
#' @examples
#' data(HMP2, package = "mgnet")
#' taxa_vars(HMP2)  
#'
#' data(subjects_HMP2, package = "mgnet")
#' taxa_vars(subjects_HMP2) 
#' 
#' @export
#' @name taxa_vars
#' @aliases taxa_vars,mgnet-method taxa_vars,mgnetList-method
setGeneric("taxa_vars", function(object) standardGeneric("taxa_vars"))

setMethod("taxa_vars", "mgnet", function(object) {
  
  if(length(object@taxa)!=0 & length(object@comm)!=0){
    return(c("comm_id", colnames(object@taxa)))
  } else if(length(object@taxa)!=0 & length(object@comm)==0){
    return(colnames(object@taxa))
  } else if(length(object@taxa)==0 & length(object@comm)!=0){
    return("comm_id")
  } else {
    return(character(length=0))
  }
  
})

setMethod("taxa_vars", "mgnetList", function(object) {
  sapply(object@mgnets, taxa_vars, simplify = FALSE, USE.NAMES = TRUE)
})

# COMM_ID
#------------------------------------------------------------------------------#
#' Retrieve Community IDs from mgnet or mgnetList Objects
#'
#' This function retrieves community IDs from an `mgnet` or `mgnetList` object, providing a flexible output
#' format. For `mgnetList`, the community IDs can be output as either a list or a structured tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt The format of the output. Possible values are:
#'   \itemize{
#'     \item `"list"`: Returns a named list where each element is a character vector of community IDs.
#'     \item `"tbl"`: Returns a tibble with columns 'mgnet', 'taxa_id', and 'comm_id', representing the 
#'           source `mgnet` object, taxa IDs, and their respective community IDs.
#'   }
#'   The default is `"list"`.
#'
#' @importFrom igraph membership
#' @importFrom tibble tibble
#' @importFrom purrr imap_dfr
#' @importFrom stats setNames
#' @export
#' @name comm_id
#' @aliases comm_id,mgnet-method comm_id,mgnetList-method
#' @examples
#' data(subject_HMP2_netw, package = "mgnet")
#' comm_id(subject_HMP2_netw)  
#'
#' data(subjects_HMP2_netw, package = "mgnet")
#' comm_id(subjects_HMP2_netw, .fmt = "list")  
#' comm_id(subjects_HMP2_netw, .fmt = "tbl")   
setGeneric("comm_id", function(object, .fmt = "list") standardGeneric("comm_id"))

setMethod("comm_id", "mgnet", function(object, .fmt = "list"){
  
  if(length(object@comm) != 0){
    
    result <- switch(.fmt,
                     list = stats::setNames(as.character(igraph::membership(object@comm)), taxa_id(object)),
                     df = data.frame("comm_id" = as.character(membership(object@comm)), row.names = taxa_id(object)),
                     tbl = tibble("taxa_id" = taxa_id(object), "comm_id" = as.character(membership(object@comm))))
    return(result)
    
  } else {
    
    return(character(0))
    
  }
})

setMethod("comm_id", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "tbl"))
  
  if (.fmt == "list") {
    
    # Return a list with each element being the community IDs from one `mgnet` object
    return(sapply(object@mgnets, function(x) comm_id(x, .fmt = "list"), 
                  simplify = FALSE, USE.NAMES = TRUE))
    
  } else {
    
    # Create a tibble with columns 'mgnet', 'taxa_id', and 'comm_id'
    results <- purrr::imap_dfr(object@mgnets, function(mgnet_obj, name) {
      tibble(
        mgnet = name,
        taxa_id = taxa_id(mgnet_obj),
        comm_id = as.character(membership(mgnet_obj@comm))
      )
    }, .id = "mgnet")
    
    return(results)
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