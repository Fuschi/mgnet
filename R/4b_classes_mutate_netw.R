#' Mutate Network Data in the mgnet Object
#'
#' This method applies transformations specified by expressions to the netw and comm slots within the mgnet object.
#' These expressions are evaluated within the context of the object or its subsets, enabling tailored manipulations based on network data.
#'
#' @param object An mgnet object.
#'        This function specifically targets the netw and comm slots, which contain network data and community structures, respectively.
#' @param ... Quoted expressions that specify how to modify the network data.
#'        The expressions can only use netw and comm as variables, which are bound to the corresponding slots in the mgnet object.
#' @param group_taxa Optional; a character vector specifying the columns to group data by before transformations.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined by group_taxa.
#'        If you do not wish to group the data leave it empty.
#'
#' ### Contextual Variables in mgnet:
#' - **netw**: Refers to the network slot in the mgnet object where network data and interactions are stored.
#' - **comm**: Refers to the community slot in the mgnet object that includes community or clustering data.
#'
#' @return Returns the modified mgnet object with updated netw and comm slots reflecting the applied transformations.
#'         The rest of the object remains unchanged, ensuring that transformations are isolated to the targeted slots.
#'
#' @export
#' @importFrom dplyr mutate select pull filter
#' @importFrom tidyr unite
#' @importFrom rlang enquos get_expr parse_expr as_quosure eval_tidy quo_is_null
#' @name mutate_netw
#' @aliases mutate_netw,mgnet-method mutate_netw,mgnetList-method
setGeneric("mutate_netw", function(object, ..., group_taxa) {standardGeneric("mutate_netw")})

setMethod("mutate_netw", "mgnet", function(object, ..., group_taxa) {
  
  if (miss_slot(object, "netw") & miss_slot(object, "comm")) {
    stop("No network and communities available in the 'mgnet' object.")
  }
  
  # Capture and check expressions
  quosures <- rlang::enquos(...)
  check_reserved_keywords(quosures)
  
  # Initialize groups with default values if it is empty.
  if(missing(group_taxa) && is_grouped_taxa(object)){
    group_taxa <- get_group_taxa(object)
  } else if(missing(group_taxa) && !is_grouped_taxa(object)){
    group_taxa <- NULL
  }
  
  # Check if is required the edge filter
  selected_links <- get_selected_links(object)
  if (!is.null(selected_links)) {
    netw0 <- netw(object)
    netw(object) <- igraph::subgraph_from_edges(graph = netw(object),
                                                eids = get_selected_links(object),
                                                delete.vertices = FALSE)
  }
  
  # Check the reserved keywords
  check_reserved_keywords(quosures)
  
  expressions_text <- purrr::map(quosures, \(x) expr_text(get_expr(x)))
  quosures_names <- names(quosures)
  quosures_names <- purrr::map2_chr(quosures_names, expressions_text,
                             \(x,y) ifelse(x == "", y, x))

  # Modify each expression
  for(i in seq_along(quosures)){
    
    text <- rlang::expr_text(rlang::get_expr(quosures[[i]]))

    if(!length(group_taxa)){
          modified_text <- str_replace_all(text, "netw", "netw(object)")
          modified_text <- str_replace_all(modified_text, "comm", "comm(object)")
    } else {
          modified_text <- str_replace_all(text, "netw", "netw(object_subset)")
          modified_text <- str_replace_all(modified_text, "comm", "comm(object_subset)")
    }
    quosures[[i]]  <- set_expr(quosures[[i]], parse_expr(modified_text))
    quosures[[i]]  <- set_env(quosures[[i]], current_env())
  }

  # Evaluate expressions and update the object
  if (!length(group_taxa)) {
    
    for (i in seq_along(quosures)) {
      taxa(object) <- gather_taxa(object) %>%
        dplyr::mutate(!!quosures_names[i] := rlang::eval_tidy(quosures[[i]], 
                                                              env = current_env()))
    }
    
  } else {
      
    subgroups <- gather_taxa(object) %>%
      dplyr::select(tidyselect::all_of(group_taxa)) %>%
      tidyr::unite("_internal_") %>%
      dplyr::pull("_internal_")
      
    unique_keys <- unique(subgroups)

    for(i in 1:length(quosures)){
      result <- vector(length = ntaxa(object))
      for(key in unique_keys){
        idx_key <- which(subgroups == key)
        object_subset <- object[, idx_key]
        result_key <- rlang::eval_tidy(quosures[[i]], env = current_env())
        result[idx_key] <- result_key
      } 
      taxa(object) <- gather_taxa(object) %>%
        mutate(!!quosures_names[i] := result)
    }
    
  }
  
  if(!is.null(selected_links)) netw(object) <- netw0
  validObject(object)
  return(object)
})


setMethod("mutate_netw", "mgnetList", function(object, ..., group_taxa) {
  
  # Capture and check expressions
  quosures <- rlang::enquos(...)
  check_reserved_keywords(quosures)
  
  # Initialize groups with default values if it is empty.
  if(missing(group_taxa) & length(get_group_taxa(object))){
    group_taxa <- get_group_taxa(object)
  } 
  
  # Check if is required the edge filter
  selected_links <- get_selected_links(object)
  if (all(all(purrr::map_lgl(selected_links, \(x) !is.null(x))))) {
    netw0 <- netw(object)
    for(mg in 1:length(object)){
      netw(object[[mg]]) <- igraph::subgraph_from_edges(graph = netw(object[[mg]]),
                                                  eids = get_selected_links(object[[mg]]),
                                                  delete.vertices = FALSE)
    }
  }
  
  # Check the reserved keywords
  check_reserved_keywords(quosures)
  
  expressions_text <- purrr::map(quosures, \(x) expr_text(get_expr(x)))
  quosures_names <- names(quosures)
  quosures_names <- purrr::map2_chr(quosures_names, expressions_text,
                             \(x,y) ifelse(x == "", y, x))
  
  # Modify each expression
  for(i in seq_along(quosures)){
    
    text <- expr_text(get_expr(quosures[[i]]))
    
    if(missing(group_taxa)){
      modified_text <- str_replace_all(text, "netw", "netw(object[[mg]])")
      modified_text <- str_replace_all(modified_text, "comm", "comm(object[[mg]])")
    } else {
      modified_text <- str_replace_all(text, "netw", "netw(object[[mg]][, idx_key])")
      modified_text <- str_replace_all(modified_text, "comm", "comm(object[[mg]][, [, idx_key]])")
    }
    quosures[[i]]  <- set_expr(quosures[[i]], parse_expr(modified_text))
    quosures[[i]]  <- set_env(quosures[[i]], current_env())
  }
  
  # Evaluate expressions and update the object
  if (missing(group_taxa)) {
    
    for (quo in seq_along(quosures)) {
      for (mg in seq_along(object)){
        taxa(object[[mg]]) <- gather_taxa(object[[mg]]) %>%
          dplyr::mutate(!!quosures_names[i] := rlang::eval_tidy(quosures[[i]], env = current_env()))
      }
    }
    
  } else {
    
    subgroups <- gather_taxa(object) %>%
      dplyr::select(tidyselect::all_of(c("mgnet", group_taxa))) %>%
      tidyr::unite("_internal_", remove = FALSE) %>%
      split(.[["mgnet"]])
    
    subgroups <- sapply(subgroups, \(x) dplyr::pull(x, "_internal_"), USE.NAMES = TRUE, simplify = FALSE)

    for(mg in seq_along(object)){
      unique_keys <- unique(subgroups[[mg]])
      for(quo in 1:length(quosures)){
        result <- vector(length = ntaxa(object[[mg]]))
        for(key in unique_keys){
          idx_key <- which(subgroups[[mg]] == key)
          result_key <- rlang::eval_tidy(quosures[[quo]], env = current_env())
          result[idx_key] <- result_key
        } 
        taxa(object[[mg]]) <- gather_taxa(object[[mg]]) %>%
          mutate(!!quosures_names[i] := result)
      }
    }
    
    
  }
  
  if (all(all(purrr::map_lgl(selected_links, \(x) !is.null(x))))) {
    for(mg in seq_along(object)) netw(object[[mg]]) <- netw0[[mg]]
  }
  validObject(object)
  return(object)
  
})



