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
#' data(mg, package = "mgnet")
#' nsample(mg)  
#'
#' data(mg, package = "mgnet")
#' nsample(mgl)  
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
#' data(mg, package = "mgnet")
#' ntaxa(mg)  
#'
#' data(mgl, package = "mgnet")
#' ntaxa(mgl)  
#'
#' @importFrom igraph vcount
#' @export
#' @name ntaxa
#' @aliases ntaxa,mgnet-method ntaxa,mgnetList-method
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))

setMethod("ntaxa", "mgnet", function(object) {
  if(length(object@abun)!=0) return(ncol(object@abun))
  else if(length(object@taxa)!=0) return(nrow(object@taxa))
  else if(length(object@rela)!=0) return(ncol(object@rela))
  else if(length(object@norm)!=0) return(ncol(object@norm))
  else if(length(object@netw)!=0) return(igraph::vcount(object@netw))
  else return(0)
})

setMethod("ntaxa", "mgnetList", function(object) {
  sapply(object@mgnets, ntaxa, simplify = TRUE, USE.NAMES = TRUE)
})


# SAMPLE_ID
#------------------------------------------------------------------------------#
#' Get sample IDs from mgnet or mgnetList Objects
#'
#'
#' @param object An `mgnet` or `mgnetList` object from which to extract sample IDs.
#' 
#' @return Retrieves the IDs of sample from an `mgnet` object or lists of sample IDs for each 
#' `mgnet` object within an `mgnetList`. 
#'         
#' @export
#' @name sample_id
#' @aliases sample_id,mgnet-method sample_id,mgnetList-method
#' 
#' @examples
#' data(mg, package = "mgnet")
#' sample_id(mg)  
#'
#' data(mgl, package = "mgnet")
#' sample_id(mgl)  
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

setMethod("sample_id", "mgnetList", function(object) {
    return(sapply(object@mgnets, sample_id, simplify = FALSE, USE.NAMES = TRUE))
})


# TAXA_ID
#------------------------------------------------------------------------------#
#' Get taxa IDs from mgnet or mgnetList Objects
#'
#'
#' @param object An `mgnet` or `mgnetList` object from which to extract taxa IDs.
#' 
#' @return Retrieves the IDs of taxa from an `mgnet` object or lists of taxa IDs for each 
#' `mgnet` object within an `mgnetList`. 
#' 
#' @importFrom igraph V
#'         
#' @export
#' @name taxa_id
#' @aliases sample_id,mgnet-method sample_id,mgnetList-method
#' 
#' @examples
#' data(mg, package = "mgnet")
#' taxa_id(mgl)  
#'
#' data(mgl, package = "mgnet")
#' taxa_id(mgl)  
#' 
#' @seealso \link{mgnet} and \link{mgnetList} for details on the classes.
setGeneric("taxa_id", function(object, .fmt = "list") standardGeneric("taxa_id"))

setMethod("taxa_id", "mgnet", function(object, .fmt = "list") {
  if(length(object@abun) != 0) return(colnames(object@abun))
  else if(length(object@taxa) != 0) return(rownames(object@taxa))
  else if(length(object@norm) != 0) return(colnames(object@norm))
  else if(length(object@rela) != 0) return(colnames(object@rela))
  else if(length(object@netw) != 0) return(igraph::V(object@netw)$name)
  else return(character(0))
})

setMethod("taxa_id", "mgnetList", function(object) {
  return(sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE))
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
#' data(mg, package = "mgnet")
#' ncomm(mg)  
#'
#' data(mgl, package = "mgnet")
#' ncomm(mgl)  
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
#' data(mg, package = "mgnet")
#' comm_id(mg)  
#'
#' data(mgl, package = "mgnet")
#' comm_id(mgl, .fmt = "list")  
#' comm_id(mgl, .fmt = "tbl")   
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


# HAS METHODS
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# HAS_SAMPLE
#------------------------------------------------------------------------------#
#' Check if Samples Are Present
#'
#' Checks whether there are samples present in an `mgnet` object or in each 
#' `mgnet` object within an `mgnetList`. For `mgnetList` objects, the output format can
#' be a list, or a single logical value indicating if any or all objects contain samples.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` contains samples, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects contain samples. Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether samples are present.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` contains samples. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects contain samples.
#'         
#' @examples
#' data(mg, package = "mgnet")
#' has_sample(mg)  
#'
#' data(mgl, package = "mgnet")
#' has_sample(mgl, .fmt = "list")  
#' has_sample(mgl, .fmt = "any")  
#' has_sample(mgl, .fmt = "all")  
#'
#' @export
#' @name has_sample
#' @aliases has_sample,mgnet-method has_sample,mgnetList-method
setGeneric("has_sample", function(object, .fmt = "list") standardGeneric("has_sample"))

setMethod("has_sample", "mgnet", function(object) {
  nsample(object) != 0
})

setMethod("has_sample", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(FALSE)
  sample_presence <- sapply(object@mgnets, has_sample, simplify = TRUE, USE.NAMES = TRUE)
  
  if (.fmt == "list") {
    return(sample_presence)
  } else if (.fmt == "any") {
    return(any(sample_presence))
  } else if (.fmt == "all") {
    return(all(sample_presence))
  }
})


# MISS_SAMPLE
#------------------------------------------------------------------------------#
#' Check if Samples are Missing in `mgnet` or `mgnetList` Objects
#'
#' This function checks whether samples are missing (i.e., absent or empty) in an `mgnet` object or 
#' in each `mgnet` object within an `mgnetList`. For `mgnetList` objects, the output format can
#' be a list, or a single logical value indicating if any or all objects are missing samples.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` is missing samples, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects are missing samples. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether samples are missing.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` is missing samples. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects are missing samples.
#'
#' @examples
#' data(mg, package = "mgnet")
#' miss_sample(mg)  
#'
#' data(mgl, package = "mgnet")
#' miss_sample(mgl, .fmt = "list")  
#' miss_sample(mgl, .fmt = "any")  
#' miss_sample(mgl, .fmt = "all")  
#'
#' @export
#' @name miss_sample
#' @aliases miss_sample,mgnet-method miss_sample-mgnetList
setGeneric("miss_sample", function(object, .fmt = "list") standardGeneric("miss_sample"))

setMethod("miss_sample", "mgnet", function(object, .fmt = "list") {
  nsample(object) == 0
})

setMethod("miss_sample", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(TRUE)
  sample_absence <- sapply(object@mgnets, miss_sample, simplify = TRUE, USE.NAMES = TRUE)
  
  if (.fmt == "list") {
    return(sample_absence)
  } else if (.fmt == "any") {
    return(any(sample_absence))
  } else if (.fmt == "all") {
    return(all(sample_absence))
  }
})


# HAS_TAXA
#------------------------------------------------------------------------------#
#' Check if Taxa Are Present
#'
#' Checks whether there are taxa present in an `mgnet` object or in each 
#' `mgnet` object within an `mgnetList`. For `mgnetList` objects, the output format can
#' be a list, or a single logical value indicating if any or all objects contain taxa.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` contains taxa, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects contain taxa. Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether taxa are present.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` contains taxa. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects contain taxa.
#'         
#' @examples
#' data(mg, package = "mgnet")
#' has_taxa(mg)  
#'
#' data(mgl, package = "mgnet")
#' has_taxa(mgl, .fmt = "list")  
#' has_taxa(mgl, .fmt = "any")  
#' has_taxa(mgl, .fmt = "all")  
#'
#' @export
#' @name has_taxa
#' @aliases has_taxa,mgnet-method has_taxa,mgnetList-method
setGeneric("has_taxa", function(object, .fmt = "list") standardGeneric("has_taxa"))

setMethod("has_taxa", "mgnet", function(object) {
  ntaxa(object) != 0
})

setMethod("has_taxa", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(FALSE)
  taxa_presence <- sapply(object@mgnets, has_taxa, simplify = TRUE, USE.NAMES = TRUE)
  
  if (.fmt == "list") {
    return(taxa_presence)
  } else if (.fmt == "any") {
    return(any(taxa_presence))
  } else if (.fmt == "all") {
    return(all(taxa_presence))
  }
})


# MISS_TAXA
#------------------------------------------------------------------------------#
#' Check if Taxa are Missing in `mgnet` or `mgnetList` Objects
#'
#' This function checks whether taxa are missing (i.e., absent or empty) in an `mgnet` object or 
#' in each `mgnet` object within an `mgnetList`. For `mgnetList` objects, the output format can
#' be a list, or a single logical value indicating if any or all objects are missing taxa.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` is missing taxa, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects are missing taxa. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether taxa are missing.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` is missing taxa. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects are missing taxa.
#'
#' @examples
#' data(mg, package = "mgnet")
#' miss_taxa(mg)  
#'
#' data(mgl, package = "mgnet")
#' miss_taxa(mgl, .fmt = "list")  
#' miss_taxa(mgl, .fmt = "any")  
#' miss_taxa(mgl, .fmt = "all")  
#'
#' @export
#' @name miss_taxa
#' @aliases miss_taxa,mgnet-method miss_taxa-mgnetList
setGeneric("miss_taxa", function(object, .fmt = "list") standardGeneric("miss_taxa"))

setMethod("miss_taxa", "mgnet", function(object, .fmt = "list") {
  length(taxa_id(object)) == 0
})

setMethod("miss_taxa", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(TRUE)
  taxa_absence <- sapply(object@mgnets, miss_taxa, simplify = TRUE, USE.NAMES = TRUE)
  
  if (.fmt == "list") {
    return(taxa_absence)
  } else if (.fmt == "any") {
    return(any(taxa_absence))
  } else if (.fmt == "all") {
    return(all(taxa_absence))
  }
})


# HAS_SLOT
#------------------------------------------------------------------------------#
#' Check if a Slot is Present in `mgnet` or `mgnetList` Objects
#'
#' This function checks whether a specific slot (e.g., `meta`, `taxa`, `netw`, etc.) 
#' is present in an `mgnet` object or in each `mgnet` object within an `mgnetList`. 
#' For `mgnetList` objects, the output format can be a list, or a single logical 
#' value indicating if any or all objects contain the specified slot.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param slot_name A character string specifying the slot to check for presence 
#'        (e.g., `"meta"`, `"taxa"`, `"netw"`, `"comm"`, `"abun"`, `"rela"`, `"norm"`).
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` contains the slot, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects contain the slot. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether the slot is present.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` contains the slot. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects contain the slot.
#'
#' @examples
#' data(mg, package = "mgnet")
#' has_slot(mg, "meta")  
#'
#' data(mgl, package = "mgnet")
#' has_slot(mgl, "taxa", .fmt = "list")  
#' has_slot(mgl, "netw", .fmt = "any")  
#' has_slot(mgl, "comm", .fmt = "all")  
#'
#' @importFrom methods slot
#' @export
#' @name has_slot
#' @aliases has_slot,mgnet-method has_slot,mgnetList-method
setGeneric("has_slot", function(object, slot_name, .fmt = "list") standardGeneric("has_slot"))

setMethod("has_slot", "mgnet", function(object, slot_name, .fmt = "list") {
  slot_value <- methods::slot(object, slot_name)
  return(length(slot_value) != 0)
})

setMethod("has_slot", "mgnetList", function(object, slot_name, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(FALSE)
  slot_presence <- sapply(object@mgnets, has_slot, slot_name = slot_name, simplify = TRUE, USE.NAMES = TRUE)
  
  if (.fmt == "list") {
    return(slot_presence)
  } else if (.fmt == "any") {
    return(any(slot_presence))
  } else if (.fmt == "all") {
    return(all(slot_presence))
  }
})


# MISS_SLOT
#------------------------------------------------------------------------------#
#' Check if a Slot is Missing in `mgnet` or `mgnetList` Objects
#'
#' This function checks whether a specific slot (e.g., `meta`, `taxa`, `netw`, etc.) 
#' is missing (i.e., empty or absent) in an `mgnet` object or in each `mgnet` object within an `mgnetList`. 
#' For `mgnetList` objects, the output format can be a list, or a single logical 
#' value indicating if any or all objects are missing the specified slot.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param slot_name A character string specifying the slot to check for absence 
#'        (e.g., `"meta"`, `"taxa"`, `"netw"`, `"comm"`, `"abun"`, `"rela"`, `"norm"`).
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` is missing the slot, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects are missing the slot. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether the slot is missing.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` is missing the slot. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects are missing the slot.
#'
#' @examples
#' data(mg, package = "mgnet")
#' miss_slot(mg, "meta")  
#'
#' data(mgl, package = "mgnet")
#' miss_slot(mgl, "taxa", .fmt = "list")  
#' miss_slot(mgl, "netw", .fmt = "any")  
#' miss_slot(mgl, "comm", .fmt = "all")  
#'
#' @export 
#' @name miss_slot
#' @aliases miss_slot,mgnet-method miss_slot,mgnetList-method
setGeneric("miss_slot", function(object, slot_name, .fmt = "list") standardGeneric("miss_slot"))

setMethod("miss_slot", "mgnet", function(object, slot_name, .fmt = "list") {
  # Dynamically check the length of the specified slot and return TRUE if missing
  slot_value <- methods::slot(object, slot_name)
  return(length(slot_value) == 0)
})

setMethod("miss_slot", "mgnetList", function(object, slot_name, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(TRUE)
  slot_absence <- sapply(object@mgnets, miss_slot, slot_name = slot_name, simplify = TRUE, USE.NAMES = TRUE)
  
  # Return the result based on the `.fmt` parameter
  if (.fmt == "list") {
    return(slot_absence)
  } else if (.fmt == "any") {
    return(any(slot_absence))
  } else if (.fmt == "all") {
    return(all(slot_absence))
  }
})


# HAS METADATA
#------------------------------------------------------------------------------#
#' Check if Taxa Metadata Are Present
#'
#' Checks whether there are taxa metadata, combining `taxa` and `comm` slots,
#' in an `mgnet` object or in each `mgnet` object within an `mgnetList`. For `mgnetList` 
#' objects, the output format can be a list, or a single logical value indicating
#' if any or all objects contain taxa.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` contains taxa, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects contain taxa. Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether taxa are present.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` contains taxa metadata. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects contain taxa.
#'         
#' @examples
#' data(mg, package = "mgnet")
#' has_metataxa(mg)  
#'
#' data(mgl, package = "mgnet")
#' has_metataxa(mgl, .fmt = "list")  
#' has_metataxa(mgl, .fmt = "any")  
#' has_metataxa(mgl, .fmt = "all")  
#'
#' @export
#' @name has_metataxa
#' @aliases has_metataxa,mgnet-method has_metataxa,mgnetList-method
setGeneric("has_metataxa", function(object, .fmt = "list") standardGeneric("has_metataxa"))

setMethod("has_metataxa", "mgnet", function(object) {
  if(length(object@taxa) != 0 || length(object@comm) != 0) TRUE else FALSE
})

setMethod("has_metataxa", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(TRUE)
  slot_presence <- sapply(object@mgnets, has_metataxa, simplify = TRUE, USE.NAMES = TRUE)
  
  # Return the result based on the `.fmt` parameter
  if (.fmt == "list") {
    return(slot_presence)
  } else if (.fmt == "any") {
    return(any(slot_presence))
  } else if (.fmt == "all") {
    return(all(slot_presence))
  }
})

# MISS METADATA
#------------------------------------------------------------------------------#
#' Check if a Taxa Metadata Missing in `mgnet` or `mgnetList` Objects
#'
#' This function checks whether taxa metadata, combining `taxa` and `comm`, 
#' is missing (i.e., empty or absent) in an `mgnet` object or in each `mgnet` object within an `mgnetList`. 
#' For `mgnetList` objects, the output format can be a list, or a single logical 
#' value indicating if any or all objects are missing the specified slot.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param slot_name A character string specifying the slot to check for absence 
#'        (e.g., `"meta"`, `"taxa"`, `"netw"`, `"comm"`, `"abun"`, `"rela"`, `"norm"`).
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of logical values, \code{"any"} for 
#'        a single logical indicating if any `mgnet` is missing the slot, and \code{"all"} 
#'        for a single logical indicating if all `mgnet` objects are missing the slot. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a logical indicating whether the slot is missing.
#'         For an `mgnetList` object, the output depends on the value of `.fmt`.
#'         If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"}, 
#'         returns a single logical value indicating if any `mgnet` is missing the slot. 
#'         If \code{.fmt = "all"}, returns a single logical value indicating if all `mgnet`
#'         objects are missing the slot.
#'
#' @examples
#' data(mg, package = "mgnet")
#' miss_metataxa(mg, "meta")  
#'
#' data(mgl, package = "mgnet")
#' miss_metataxa(mgl, "taxa", .fmt = "list")  
#' miss_metataxa(mgl, "netw", .fmt = "any")  
#' miss_metataxa(mgl, "comm", .fmt = "all")  
#'
#' @export 
#' @name miss_metataxa
#' @aliases miss_metataxa,mgnet-method miss_metataxa,mgnetList-method
setGeneric("miss_metataxa", function(object, .fmt = "list") standardGeneric("miss_metataxa"))

setMethod("miss_metataxa", "mgnet", function(object) {
  if(length(object@taxa) == 0 && length(object@comm) == 0) TRUE else FALSE
})

setMethod("miss_metataxa", "mgnetList", function(object, .fmt = "list") {
  .fmt <- match.arg(.fmt, choices = c("list", "any", "all"))
  
  if(length(object) == 0) return(TRUE)
  slot_absence <- sapply(object@mgnets, miss_metataxa, simplify = TRUE, USE.NAMES = TRUE)
  
  # Return the result based on the `.fmt` parameter
  if (.fmt == "list") {
    return(slot_absence)
  } else if (.fmt == "any") {
    return(any(slot_absence))
  } else if (.fmt == "all") {
    return(all(slot_absence))
  }
})