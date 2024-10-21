#' Rank Taxa Based on Custom Criteria in `mgnet` Objects
#'
#' This function sorts and ranks taxa based on a specified aggregation or calculation 
#' applied to abundance-related variables (`abun`, `rela`, `norm`) or other grouping 
#' functions. It is designed for `mgnet` or `mgnetList` objects to facilitate focused 
#' analyses on grouped or aggregated data.
#'
#' @param object An `mgnet` object.
#' @param expression A dplyr-compatible expression that defines the criteria for 
#'        sorting taxa. This should involve only the abundance variables `abun`, `rela`, 
#'        `norm` or grouping functions from `dplyr`. The expression should not reference 
#'        other data keys directly.
#' @param by Optional character string specifying the taxa metadata column to group by 
#'        before applying the ranking expression. Default is `"taxa_id"`.
#' @param decreasing Logical indicating whether sorting should be in decreasing order. 
#'        Default is `TRUE`.
#' @param mode For `mgnetList` objects, specifies whether to evaluate the expression on the merged objects (`"merged"`) 
#'        or evaluate each object separately (`"separate"`). This parameter is ignored for `mgnet` objects.
#'
#' @return A vector containing the `by` column, sorted according to the specified 
#'         `expression`.
#'
#' @details
#' `top_taxa` allows users to specify complex sorting and ranking expressions that 
#' are applied within groups defined by the `by` parameter. The function checks if the 
#' provided expression uses only allowed variables (`abun`, `rela`, `norm`) and grouping 
#' functions like `n()`, ensuring the expression's appropriateness for the operation.
#'
#' This function dynamically constructs a data frame from the `mgnet` object's slots,
#' performs the necessary transformations and aggregations, and then sorts the results.
#' If no `abun`, `rela`, or `norm` keys are needed for the expression, it simply sorts the
#' taxa metadata. Otherwise, it performs left joins to include the necessary abundance data.
#'
#' @importFrom dplyr group_by left_join summarise pull
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom methods slot
#' @export
#' @aliases top_taxa,mgnet-method top_taxa,mgnetList-method
setGeneric("top_taxa", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged") standardGeneric("top_taxa"))

setMethod("top_taxa", "mgnet", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged"){
  
  # Capture the expression and ensure it is a single, simple expression
  expr_quo <- enquo(expression)
  
  # Extract variable names from the quosure expression
  keys_required <- all.vars(expr_quo)
  
  # # Check keys properties
  # if(length(keys_required) > 0){
  #   if(!all(keys_required %in% c("abun", "rela", "norm"))){
  #     stop("Error: The expression for sorting taxa should only involve 'abun', 'rela', 'norm' or grouping functions from `dplyr` such as 'n()' without referencing other data keys.")
  #   }
  # }
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun", "rela", "norm"))
  
  # Validate the 'by' argument
  if(is.null(by)) {
    by <- "taxa_id"
  }
  
  available_taxa_vars <- c("taxa_id", taxa_vars(object))
  if (!is.character(by) || length(by) != 1 || !by %in% available_taxa_vars) {
    stop(paste0("Error: 'by' must be a single string in this list: ", toString(available_taxa_vars)))
  }
  
  # # Forbidden functions and disallowed variables
  check_forbidden_expressions(enquos(expression))

  # END CHECKS
  #----------------------------------------------------------------------------#

  # CONSTRUCT THE BASE DATA.FRAME FOR THE EVALUATION
  #----------------------------------------------------------------------------#
  if(length(needed_abundance_keys)==0){
    
    if (length(taxa(object)) != 0) {
      taxa_sorted <- taxa(object, .fmt = "tbl") 
    } else {
      taxa_sorted <- tibble::tibble(taxa_id = taxa_id(object))
    }
    
  } else {
    
    taxa_sorted <- tidyr::expand_grid(sample_id = sample_id(object), taxa_id = taxa_id(object))
    
    # Left join all abundance variables needed
    for (abundance_key in needed_abundance_keys) {
      taxa_sorted <- taxa_sorted %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id,
                                names_to = "taxa_id",
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
    
    taxa_sorted <- taxa_sorted %>%
      dplyr::left_join(taxa(object, "tbl"), by = "taxa_id")
    
  }
  
  taxa_sorted <- taxa_sorted %>%
    dplyr::group_by(!!!rlang::syms(by)) %>%
    dplyr::summarise(!!expr_quo, .groups = "drop") %>%
    as.data.frame()

  taxa_sorted <- taxa_sorted[order(taxa_sorted[, ncol(taxa_sorted)], decreasing = decreasing), ]
  
  
  return(pull(taxa_sorted, {{by}}))
})

setMethod("top_taxa", "mgnetList", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged"){
  
  mode <- match.arg(mode, choices = c("merged", "separate"))
  
  # Capture the expression and ensure it is a single, simple expression
  expr_quo <- enquo(expression)
  
  # Extract variable names from the quosure expression
  keys_required <- all.vars(expr_quo)
  
  # Check keys properties
  if(length(keys_required) > 0){
    if(!all(keys_required %in% c("abun", "rela", "norm"))){
      stop("Error: The expression for sorting taxa should only involve 'abun', 'rela', 'norm' or grouping functions from `dplyr` such as 'n()' without referencing other data keys.")
    }
  }
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun", "rela", "norm"))
  
  # Validate the 'by' argument
  if(is.null(by)) {
    by <- "taxa_id"
  }
  
  # # Forbidden functions and disallowed variables
  check_forbidden_expressions(enquos(expression))
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # IMPLEMENTATION
  #----------------------------------------------------------------------------#
  if(mode == "merged"){
    
    # MERGED
    #--------------------------------------------------------#
    
    available_taxa_vars <- c("taxa_id", unique(unlist(taxa_vars(object))) )
    if (!is.character(by) || length(by) != 1 || !by %in% available_taxa_vars) {
      stop(paste0("Error: 'by' must be a single string in this list: ", toString(available_taxa_vars)))
    }
    
    if(length(needed_abundance_keys)==0){
      
      taxa_sorted_merged <- object %>%
        purrr::map(\(x) {
          if(length(taxa(x)) != 0){
            taxa(x, .fmt = "tbl")
          } else {
            tibble::tibble(taxa_id = taxa_id(x))
          }
        }) %>%
        purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
        purrr::list_rbind()
      
    } else {
      
      if(length(needed_abundance_keys)!=0){
        
        taxa_sorted_merged <- purrr::map(object, \(mgnet_obj){
          
          taxa_sorted <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                            taxa_id = taxa_id(mgnet_obj))
          
          for(abundance_key in needed_abundance_keys){
            taxa_sorted <- taxa_sorted %>%
              dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                                 tibble::as_tibble(rownames = "sample_id") %>%
                                 tidyr::pivot_longer(-sample_id, 
                                                     names_to = "taxa_id", 
                                                     values_to = abundance_key),
                               by = c("sample_id", "taxa_id"))
          }
          return(taxa_sorted)
        }) %>%
          purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
          purrr::list_rbind()
        
        if(nrow(taxa(object, "tbl")) > 0){
          taxa_sorted_merged <- taxa_sorted_merged %>%
            left_join(taxa(object, "tbl"), dplyr::join_by(mgnet, taxa_id))
        }
      }
    }
    
    taxa_sorted_merged <- taxa_sorted_merged %>%
      group_by(!!!rlang::syms(by)) %>%
      reframe(!!expr_quo) %>%
      as.data.frame()
    
    taxa_sorted_merged <- taxa_sorted_merged[order(taxa_sorted_merged[, ncol(taxa_sorted_merged)], decreasing = decreasing), ]
    
    
    return(dplyr::pull(taxa_sorted_merged, {{by}}))
    
  } else if(mode == "separate"){
    
    # Check that 'by' is present in all mgnets
    for(mgnet in object@mgnets) {
      if(!(by %in% colnames(taxa(mgnet, "tbl")))) {
        stop(paste("Error: 'by' variable", by, "is not present in all mgnet objects' taxa metadata."))
      }
    }
    
    # SEPARATE
    #--------------------------------------------------------#
    separate_sorted <- sapply(object, \(x){
      
      if(length(needed_abundance_keys)==0){
        
        taxa_sorted <- if(length(taxa(x))!=0){
          taxa(x, .fmt = "tbl")
        } else {
          tibble::tibble(taxa_id = taxa_id(x))
        }
        
      } else {
        
        taxa_sorted <- tidyr::expand_grid(sample_id = sample_id(x),
                                          taxa_id = taxa_id(x))
        
        for(abundance_key in needed_abundance_keys){
          taxa_sorted <- taxa_sorted %>%
            dplyr::left_join(methods::slot(x, abundance_key) %>%
                               tibble::as_tibble(rownames = "sample_id") %>%
                               tidyr::pivot_longer(-sample_id, 
                                                   names_to = "taxa_id", 
                                                   values_to = abundance_key),
                             by = c("sample_id", "taxa_id"))
        }
          
        if(nrow(taxa(x, "tbl")) > 0){
          taxa_sorted <- taxa_sorted %>%
            left_join(taxa(x, "tbl"), by = "taxa_id")
        }
      }
      
      taxa_sorted <- taxa_sorted %>%
        group_by(!!!rlang::syms(by)) %>%
        dplyr::summarise(!!expr_quo, .groups = "drop") %>%
        as.data.frame()
      
      taxa_sorted <- taxa_sorted[order(taxa_sorted[, ncol(taxa_sorted)], decreasing = decreasing), ]
      
      
      return(dplyr::pull(taxa_sorted, {{by}}))
      
    }, simplify = FALSE, USE.NAMES = TRUE)
    return(separate_sorted)
  }
  
})


#' Rank Sample Based on Custom Criteria in `mgnet` Objects
#'
#' This function sorts and ranks taxa based on a specified aggregation or calculation 
#' applied to abundance-related variables (`abun`, `rela`, `norm`) or other grouping 
#' functions. It is designed for `mgnet` or `mgnetList` objects to facilitate focused 
#' analyses on grouped or aggregated data.
#'
#' @param object An `mgnet` object.
#' @param expression A dplyr-compatible expression that defines the criteria for 
#'        sorting samples. This should involve only the abundance variables `abun`, `rela`, 
#'        `norm` or grouping functions from `dplyr`. The expression should not reference 
#'        other data keys directly.
#' @param by Optional character string specifying the sample metadata column to group by 
#'        before applying the ranking expression. Default is `"sample_id"`.
#' @param decreasing Logical indicating whether sorting should be in decreasing order. 
#'        Default is `TRUE`.
#' @param mode For `mgnetList` objects, specifies whether to evaluate the expression on the merged objects (`"merged"`) 
#'        or evaluate each object separately (`"separate"`). This parameter is ignored for `mgnet` objects.
#'
#' @return A vector containing the `by` column, sorted according to the specified 
#'         `expression`.
#'
#' @details
#' `top_sample` allows users to specify complex sorting and ranking expressions that 
#' are applied within groups defined by the `by` parameter. The function checks if the 
#' provided expression uses only allowed variables (`abun`, `rela`, `norm`) and grouping 
#' functions like `n()`, ensuring the expression's appropriateness for the operation.
#'
#' This function dynamically constructs a data frame from the `mgnet` object's slots,
#' performs the necessary transformations and aggregations, and then sorts the results.
#' If no `abun`, `rela`, or `norm` keys are needed for the expression, it simply sorts the
#' taxa metadata. Otherwise, it performs left joins to include the necessary abundance data.
#'
#' @importFrom dplyr group_by left_join summarise pull
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom methods slot
#' @export
#' @aliases top_sample,mgnet-method top_sample,mgnetList-method
setGeneric("top_sample", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged") standardGeneric("top_sample"))

setMethod("top_sample", "mgnet", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged"){
  
  # Capture the expression and ensure it is a single, simple expression
  expr_quo <- enquo(expression)
  
  # Extract variable names from the quosure expression
  keys_required <- all.vars(expr_quo)
  
  # Check keys properties
  if(length(keys_required) > 0){
    if(!all(keys_required %in% c("abun", "rela", "norm"))){
      stop("Error: The expression for sorting taxa should only involve 'abun', 'rela', 'norm' or grouping functions from `dplyr` such as 'n()' without referencing other data keys.")
    }
  }
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun", "rela", "norm"))
  
  # Validate the 'by' argument
  if(is.null(by)) {
    by <- "sample_id"
  }
  
  available_meta_vars <- c("sample_id", meta_vars(object))
  if (!is.character(by) || length(by) != 1 || !by %in% available_meta_vars) {
    stop(paste0("Error: 'by' must be a single string in this list: ", toString(available_meta_vars)))
  }
  
  # # Forbidden functions and disallowed variables
  check_forbidden_expressions(enquos(expression))
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CONSTRUCT THE BASE DATA.FRAME FOR THE EVALUATION
  #----------------------------------------------------------------------------#
  if(length(needed_abundance_keys)==0){
    
    if (length(meta(object)) != 0) {
      meta_sorted <- meta(object, .fmt = "tbl") 
    } else {
      meta_sorted <- tibble::tibble(sample_id = sample_id(object))
    }
    
  } else {
    
    meta_sorted <- tidyr::expand_grid(sample_id = sample_id(object), taxa_id = taxa_id(object))
    
    # Left join all abundance variables needed
    for (abundance_key in needed_abundance_keys) {
      meta_sorted <- meta_sorted %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id,
                                names_to = "taxa_id",
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
    
    meta_sorted <- meta_sorted %>%
      dplyr::left_join(meta(object, "tbl"), by = "sample_id")
    
  }
  
  meta_sorted <- meta_sorted %>%
    dplyr::group_by(!!!rlang::syms(by)) %>%
    dplyr::summarise(!!expr_quo, .groups = "drop") %>%
    as.data.frame()
  
  meta_sorted <- meta_sorted[order(meta_sorted[, ncol(meta_sorted)], decreasing = decreasing), ]
  
  
  return(pull(meta_sorted, {{by}}))
})

setMethod("top_sample", "mgnetList", function(object, expression, by = NULL, decreasing = TRUE, mode = "merged"){
  
  mode <- match.arg(mode, choices = c("merged", "separate"))
  
  # Capture the expression and ensure it is a single, simple expression
  expr_quo <- enquo(expression)
  
  # Extract variable names from the quosure expression
  keys_required <- all.vars(expr_quo)
  
  # Check keys properties
  if(length(keys_required) > 0){
    if(!all(keys_required %in% c("abun", "rela", "norm"))){
      stop("Error: The expression for sorting taxa should only involve 'abun', 'rela', 'norm' or grouping functions from `dplyr` such as 'n()' without referencing other data keys.")
    }
  }
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun", "rela", "norm"))
  
  # Validate the 'by' argument
  if(is.null(by)) {
    by <- "sample_id"
  }
  
  # # Forbidden functions and disallowed variables
  check_forbidden_expressions(enquos(expression))
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # IMPLEMENTATION
  #----------------------------------------------------------------------------#
  if(mode == "merged"){
    
    available_meta_vars <- c("sample_id", unique(unlist(meta_vars(object))) )
    if (!is.character(by) || length(by) != 1 || !by %in% available_meta_vars) {
      stop(paste0("Error: 'by' must be a single string in this list: ", toString(available_meta_vars)))
    }
    
    # UNITE
    #--------------------------------------------------------#
    if(length(needed_abundance_keys)==0){
      
      meta_sorted_merged <- object %>%
        purrr::map(\(x) {
          if(length(meta(x)) != 0){
            meta(x, .fmt = "tbl")
          } else {
            tibble::tibble(taxa_id = taxa_id(x))
          }
        }) %>%
        purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
        purrr::list_rbind()
      
    } else {
      
      if(length(needed_abundance_keys)!=0){
        
        meta_sorted_merged <- purrr::map(object, \(mgnet_obj){
          
          meta_sorted <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                            taxa_id = taxa_id(mgnet_obj))
          
          for(abundance_key in needed_abundance_keys){
            meta_sorted <- meta_sorted %>%
              dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                                 tibble::as_tibble(rownames = "sample_id") %>%
                                 tidyr::pivot_longer(-sample_id, 
                                                     names_to = "taxa_id", 
                                                     values_to = abundance_key),
                               by = c("sample_id", "taxa_id"))
          }
          return(meta_sorted)
        }) %>%
          purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
          purrr::list_rbind()
        
        if(nrow(meta(object, "tbl")) > 0){
          meta_sorted_merged <- meta_sorted_merged %>%
            left_join(meta(object, "tbl"), dplyr::join_by(mgnet, sample_id))
        }
      }
    }
    
    meta_sorted_merged <- meta_sorted_merged %>%
      group_by(!!!rlang::syms(by)) %>%
      reframe(!!expr_quo) %>%
      as.data.frame()
    
    meta_sorted_merged <- meta_sorted_merged[order(meta_sorted_merged[, ncol(meta_sorted_merged)], decreasing = decreasing), ]
    
    
    return(dplyr::pull(meta_sorted_merged, {{by}}))
    
  } else if(mode == "separate"){
    
    # Check that 'by' is present in all mgnets
    for(mgnet in object@mgnets) {
      if(!(by %in% colnames(meta(mgnet, "tbl")))) {
        stop(paste("Error: 'by' variable", by, "is not present in all mgnet objects' sample metadata."))
      }
    }
    
    # SEPARATE
    #--------------------------------------------------------#
    separate_sorted <- sapply(object, \(x){
      
      if(length(needed_abundance_keys)==0){
        
        meta_sorted <- if(length(meta(x))!=0){
          meta(x, .fmt = "tbl")
        } else {
          tibble::tibble(sample_id = sample_id(x))
        }
        
      } else {
        
        meta_sorted <- tidyr::expand_grid(sample_id = sample_id(x),
                                          taxa_id = taxa_id(x))
        
        for(abundance_key in needed_abundance_keys){
          meta_sorted <- meta_sorted %>%
            dplyr::left_join(methods::slot(x, abundance_key) %>%
                               tibble::as_tibble(rownames = "sample_id") %>%
                               tidyr::pivot_longer(-sample_id, 
                                                   names_to = "taxa_id", 
                                                   values_to = abundance_key),
                             by = c("sample_id", "taxa_id"))
        }
        
        if(nrow(meta(x, "tbl")) > 0){
          meta_sorted <- meta_sorted %>%
            left_join(meta(x, "tbl"), by = "sample_id")
        }
      }
      
      meta_sorted <- meta_sorted %>%
        group_by(!!!rlang::syms(by)) %>%
        dplyr::summarise(!!expr_quo, .groups = "drop") %>%
        as.data.frame()
      
      meta_sorted <- meta_sorted[order(meta_sorted[, ncol(meta_sorted)], decreasing = decreasing), ]
      
      
      return(dplyr::pull(meta_sorted, {{by}}))
      
    }, simplify = FALSE, USE.NAMES = TRUE)
    return(separate_sorted)
  }
  
})