#------------------------------------------------------------------------------#
#' Apply Mutations with Abundance Check
#'
#' This internal function applies mutations to a data frame, handling both 
#' abundance-related and non-abundance-related expressions.
#'
#' @param info_mutated A data frame representing the metadata.
#' @param long_abun A data frame containing abundance-related data, 
#'        used for abundance-related expressions.
#' @param expressions A list of expressions to apply to the data.
#' @param .by A character vector indicating the grouping variable(s).
#' @param sample_or_taxa ...
#'
#' @return A mutated version of the input `taxa_mutated_merged` data frame.
#' @keywords internal
apply_mutate_verb <- function(info_mutated, long_abun, expressions, .by, sample_or_taxa) {
  
  join_by <- if(sample_or_taxa=="sample") "sample_id" else "taxa_id"
  del_key <- if(sample_or_taxa=="sample") "taxa_id" else "sample_id"
  if("mgnet" %in% colnames(info_mutated)) join_by <- c("mgnet", join_by)
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      info_mutated <- long_abun %>%
        dplyr::left_join(info_mutated, by = join_by) %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c({{del_key}}, "abun", "rela", "norm"))) %>%
        dplyr::distinct() 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      info_mutated <- info_mutated %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  } # loop expressions
  
  return(info_mutated)
}

#------------------------------------------------------------------------------#
#' Apply Filters with Abundance Check
#'
#' This internal function applies filters to a data frame, handling both 
#' abundance-related and non-abundance-related expressions.
#'
#' @param info_mutated A data frame representing the metadata.
#' @param long_abun A data frame containing abundance-related data, 
#'        used for abundance-related expressions.
#' @param expressions A list of expressions to apply to the data.
#' @param .by A character vector indicating the grouping variable(s).
#' @param sample_or_taxa ...
#'
#' @return A mutated version of the input `info_mutated` data frame.
#' @keywords internal
apply_filter_verb <- function(metadata, long_abun, expressions, .by, sample_or_taxa) {
  
  join_by <- if(sample_or_taxa=="sample") "sample_id" else "taxa_id"
  del_key <- if(sample_or_taxa=="sample") "taxa_id" else "sample_id"
  pull_id <- if(sample_or_taxa=="sample") "sample_id" else "taxa_id"
  
  if("mgnet" %in% colnames(metadata)) join_by <- c("mgnet", join_by)
  filtered_id <- list()
  
  for (i in seq_along(expressions)) {
    
    expr_vars <- all.vars(expressions[[i]])
    
    if (any(expr_vars %in% c("abun","rela","norm"))) {
      
      # Process expressions involving abundance-related variables
      filtered_id[[i]] <- long_abun %>%
        dplyr::left_join(metadata, by = {{join_by}}) %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c({{del_key}}, "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::pull({{pull_id}})
      
    } else {
      
      # Process expressions not involving abundance-related variables
      filtered_id[[i]] <- metadata %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::filter(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::pull({{pull_id}})
      
    }
  }
  
  return(filtered_id)
  
}


#------------------------------------------------------------------------------#
#' @keywords internal
apply_expr_var_sort <- function(metadata, long_abun, expression, var, n, decreasing) {
      
    expr_vars <- all.vars(expression)

    if (any(expr_vars %in% c("abun", "rela", "norm"))) {

      var_sorted <- long_abun %>%
        dplyr::group_by(!!!rlang::syms(var)) %>%
        dplyr::reframe('_internal_' := !!expression) %>%
        as.data.frame()
      
    } else {

      var_sorted <- metadata %>%
        dplyr::group_by(!!!rlang::syms(var)) %>%
        dplyr::reframe('_internal_' := !!expression) %>%
        as.data.frame()

    }
    
    if(!is.logical(var_sorted[, '_internal_']) && !is.numeric(var_sorted[, '_internal_']) || !is.vector(var_sorted[, '_internal_'])){
      stop("Error: The results of sorting expressions must be numeric or logical vector.")
    }
    
    if(is.logical(var_sorted[, '_internal_'])){
      n <- sum(var_sorted[, '_internal_'])
      var_sorted <- var_sorted[order(var_sorted[, '_internal_'], decreasing = decreasing), ]
    } else {
      n <- min(n, nrow(var_sorted))
      var_sorted <- var_sorted[order(var_sorted[, '_internal_'], decreasing = TRUE), ]
    }

    return(list(var_sorted = var_sorted, n = n))
}







