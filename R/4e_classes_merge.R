#' Merge Multiple `mgnet` or `mgnetList` Objects
#'
#' @description
#' Merges multiple `mgnet` or `mgnetList` objects into a single `mgnet` or `mgnetList` object by merging their respective slots according to user-specified functions. This method provides flexibility to define how each slot is combined, which is essential because a default merging behavior for such complex data structures is not feasible.
#'
#' @param ... `mgnet` objects, `mgnetList` objects, or lists of such objects to be merged.
#' @param abun A function to merge the 'abun' (abundance) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param rela A function to merge the 'rela' (relative abundance) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param norm A function to merge the 'norm' (normalized abundance) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param meta A function to merge the 'meta' (metadata) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param taxa A function to merge the 'taxa' (taxonomic data) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param netw A function to merge the 'netw' (network) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#' @param comm A function to merge the 'comm' (community structure) slots from the provided `mgnet` objects, or `NULL` if no merging is required for this slot.
#'
#' @return A new `mgnet` or `mgnetList` object containing the merged data from the input objects. Each slot of the returned object is the result of the merging function applied to that slot, or remains `NULL` if no function was specified for that slot.
#'
#' @details
#' The function allows specifying custom merging strategies for each slot:
#' - **Abundance data (`abun`)**: Functions like `dplyr::bind_cols` or `dplyr::bind_rows` might be appropriate.
#' - **Metadata (`meta`)**: Merging can be handled by functions such as `dplyr::full_join` or `dplyr::bind_rows`.
#' - **Network structures (`netw`)**: These might be combined using `igraph::union` or similar graph-specific functions.
#' This approach ensures flexibility in handling the complex and varied data stored within `mgnet` objects.
#'
#' The method checks for consistency among `mgnetList` objects. If merging `mgnetList` objects, it first ensures that all have identical names in the 'mgnets' slot. Mismatched names will cause the function to stop and issue an error.
#'
#' @importFrom methods new
#' @importFrom purrr map reduce map_lgl list_flatten
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#'   # Assuming mg1 and mg2 are mgnet objects
#'   merged_mgnet <- merge_mgnets(mg1, mg2, abun = rbind)
#'   # For merging mgnetList objects with custom functions for each slot
#'   merged_list <- merge_mgnets(list1, list2, 
#'                               meta = function(x) dplyr::bind_rows(x),
#'                               netw = igraph::graph.union)
#' }
merge_mgnets <- function(..., 
                         abun = NULL, rela = NULL, norm = NULL,
                         meta = NULL, taxa = NULL,
                         netw = NULL, comm = NULL) {
  
  # Flatten input to handle lists of mgnet objects or mgnetList
  mgnets <- purrr::list_flatten(list(...))
  
  if(length(mgnets) <= 1) return(mgnets)
  
  # Check for non-mgnet objects and stop execution if found
  are_mgnet     <- all(purrr::map_lgl(mgnets, is_mgnet))
  are_mgnetList <- all(purrr::map_lgl(mgnets, is_mgnetList))
  if (isFALSE(are_mgnet) & isFALSE(are_mgnetList)) {
    stop("Error: All elements must be 'mgnet' or 'mgnetList'.")
  }
  
  # Define a helper function to check if all items are equal
  all_equal <- function(items) {
    all(sapply(2:length(items), function(i) identical(items[[1]], items[[i]])))
  }
  
  # MGNET
  #--------------------------------------------#
  if(are_mgnet){
    new_mgnet <- new("mgnet")
    
    # Merge each slot using the provided functions or leave them as NULL
    if(!is.null(abun)){
      new_mgnet@abun <- purrr::map(mgnets, \(x) mgnet::abun(x, .fmt = "df")) %>%
        purrr::reduce(abun) %>%
        as.matrix()
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::abun(x)))){
      new_mgnet@abun <- mgnet::abun(mgnets[[1]])
    }
    
    if(!is.null(rela)){
      new_mgnet@rela <- purrr::map(mgnets, \(x) mgnet::rela(x, .fmt = "df")) %>%
        purrr::reduce(rela) %>%
        as.matrix()
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::rela(x)))){
      new_mgnet@rela <- mgnet::rela(mgnets[[1]])
    }
    
    if(!is.null(norm)){
      new_mgnet@norm <- purrr::map(mgnets, \(x) mgnet::norm(x, .fmt = "df")) %>%
        purrr::reduce(norm) %>%
        as.matrix()
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::norm(x)))){
      new_mgnet@norm <- mgnet::norm(mgnets[[1]])
    }
    
    if(!is.null(meta)){
      new_mgnet@meta <- purrr::map(mgnets, \(x) mgnet::meta(x)) %>%
        purrr::reduce(meta)
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::meta(x)))){
      new_mgnet@meta <- mgnet::meta(mgnets[[1]])
    }
    
    if(!is.null(taxa)){
      new_mgnet@taxa <- purrr::map(mgnets, \(x) mgnet::taxa(x)) %>%
        purrr::reduce(taxa)
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::taxa(x)))){
      new_mgnet@taxa <- mgnet::taxa(mgnets[[1]])
    }
    
    if(!is.null(netw)){
      new_mgnet@netw <- purrr::map(mgnets, \(x) mgnet::netw(x)) %>%
        purrr::reduce(netw)
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::netw(x)))){
      new_mgnet@netw <- mgnet::netw(mgnets[[1]])
    }
    
    if(!is.null(comm)){
      new_mgnet@comm <- purrr::map(mgnets, \(x) mgnet::comm(x)) %>%
        purrr::reduce(comm)
    } else if(all_equal(purrr::map(mgnets, \(x) mgnet::comm(x)))){
      new_mgnet@comm <- mgnet::comm(mgnets[[1]])
    }
    
    validObject(new_mgnet)
    return(new_mgnet)
  }
  
  
  # MGNETLIST
  #--------------------------------------------#
  if(are_mgnetList){
    
    names_list <- purrr::map(mgnets, names)
    
    # Check if all lists of names are equal
    are_names_equal <- purrr::reduce(names_list, \(a,b) all(a%in%b) & all(b%in%a))
    
    # If names are not equal, stop the execution and report the inconsistency
    if (!are_names_equal) {
      stop("Error: The mgnetList objects have different names in the 'mgnets' slot.")
    }
    
    unique_names <- unique(unlist(names_list))
    merged_list <- purrr::map(unique_names, \(name) {
      
      mgnets_name <- purrr::map(mgnets, \(x){x[[name]]})
      return(merge_mgnets(mgnets_name,
                          abun = abun, rela = rela, norm = norm,
                          meta = meta, taxa = taxa,
                          netw = netw, comm = comm))
    }) %>% 
      stats::setNames(unique_names) %>% 
      mgnetList()
    return(merged_list)
  }
  
}


#' Collapse Specified Combinations within an mgnetList
#'
#' @description
#' Merges specific combinations of `mgnet` objects within an `mgnetList` according to user-defined rules,
#' allowing for targeted merging of subsets within the list. This function facilitates focused analyses
#' by enabling the aggregation of related data into meaningful groups.
#'
#' @param object An `mgnetList` containing multiple `mgnet` objects.
#' @param by A list specifying the combinations of `mgnet` object names to be merged. Each element in
#'           the list should be a character vector containing the names of the `mgnet` objects to merge.
#'           If NULL, all `mgnet` objects in the list are merged into a single `mgnet` object.
#' @param abun A function to merge the 'abun' (abundance) slots from the selected `mgnet` objects, or
#'             NULL if no merging is required for this slot.
#' @param rela A function to merge the 'rela' (relative abundance) slots, or NULL if not required.
#' @param norm A function to merge the 'norm' (normalized abundance) slots, or NULL if not required.
#' @param meta A function to merge the 'meta' (metadata) slots, or NULL if not required.
#' @param taxa A function to merge the 'taxa' (taxonomic data) slots, or NULL if not required.
#' @param netw A function to merge the 'netw' (network) slots, or NULL if not required.
#' @param comm A function to merge the 'comm' (community structure) slots, or NULL if not required.
#'
#' @return An `mgnetList` containing the merged `mgnet` objects. Each entry in the returned `mgnetList`
#'         corresponds to a merged group as defined in the `by` parameter. The names of the entries
#'         reflect the groups specified in `by` or are derived by concatenating the names of the merged
#'         `mgnet` objects with dashes if `by` is unnamed.
#'
#' @details
#' The function checks for the validity of names specified in the `by` parameter against the names
#' of `mgnet` objects in the `mgnetList`. It ensures that all specified groups contain valid `mgnet`
#' object names and that each group's name is unique. This function is useful for scenarios where
#' specific subsets of data within an `mgnetList` need to be aggregated based on certain criteria
#' or experimental conditions.
#'
#' If `by` is NULL, a default merge of all `mgnet` objects in the list is performed, resulting in a
#' single `mgnet` object.
#'
#' Example usage:
#' \dontrun{
#'   # Assuming list1 contains several mgnet objects named 'A', 'B', 'C', etc.
#'   result <- mgnet_collapse(list1, by = list(Group1 = c("A", "B"), Group2 = c("C")),
#'                            meta = function(x) dplyr::bind_rows(x),
#'                            netw = igraph::graph.union)
#' }
#'
#' @importFrom stats setNames
#' @name mgnet_collapse
#' @aliases mgnet_collapse,mgnetList-method
#' @export
setGeneric("mgnet_collapse", function(object, by = NULL,
                                      abun = NULL, rela = NULL, norm = NULL, 
                                      meta = NULL, taxa = NULL, 
                                      netw = NULL, comm = NULL) {
  standardGeneric("mgnet_collapse")
})

setMethod("mgnet_collapse", "mgnetList", function(object, by = NULL,
                                                  abun = NULL, rela = NULL, norm = NULL, 
                                                  meta = NULL, taxa = NULL, 
                                                  netw = NULL, comm = NULL) {
  
  if(is.null(by)){
    return(merge_mgnets(as.list(object), 
                        abun = abun, rela = rela, norm = norm,
                        meta = meta, taxa = taxa, 
                        netw = netw, comm = comm))
  } 
  
  if (!is.list(by)) {
    stop("'by' must be a list of character vectors specifying combinations of mgnet names.")
  }
  
  # Check the validity of elements in 'by'
  if (!all(sapply(by, \(x) is.character(x) && all(x %in% names(object))))) {
    stop("All elements in 'by' must be character vectors containing valid names of mgnet objects within the mgnetList.")
  }
  
  # Prepare names for the output based on 'by' or create descriptive names
  names(by) <- ifelse(nzchar(names(by)), names(by), sapply(by, \(x) paste(x, collapse = "-")))
  
  # Merge specified mgnet combinations
  results <- stats::setNames(vector("list", length(by)), names(by))
  for (i in seq_along(by)) {
    sub_object <- mgnets(object)[by[[i]]]
    results[[i]] <- merge_mgnets(sub_object, 
                                 abun = abun, rela = rela, norm = norm,
                                 meta = meta, taxa = taxa, 
                                 netw = netw, comm = comm)
  }
  
  return(mgnetList(results))
})

