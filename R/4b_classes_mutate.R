#' Modify and Augment `mgnet` Objects by Transforming the `meta` Slot
#'
#' This function dynamically manipulates the `sample` slot within `mgnet` or `mgnetList` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        The function targets the `sample` slot, which contains metadata for each sample.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the sample slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'sample_id' for `mgnet` objects and c('mgnet', 'sample_id') for `mgnetList` objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by `.by`. Usage of 'taxa_id' as a grouping variable is strictly prohibited to maintain
#'        data consistency and avoid misinterpretation.
#'
#' @details The function is designed to integrate seamlessly with the `tidyverse`, allowing users
#'          to utilize familiar and potent data manipulation verbs such as `mutate`, `filter`.
#'          It supports using any `tidyverse`-compatible expressions, including conditional operations,
#'          summarizations, and complex transformations involving both abundance and metadata fields.
#'          This flexibility makes it particularly useful for ecological and biological data analysis,
#'          where combining different data types and conditions is common.
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **sample_id**: An essential identifier used to uniquely reference individual samples within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnetList` object with updated `sample` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_meta,mgnet-method mutate_meta,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate arrange
#' @importFrom tidyr expand_grid 
#' @importFrom tidyselect any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_meta", function(object, ..., .by) {standardGeneric("mutate_meta")})

setMethod("mutate_meta", "mgnet", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (nsample(object) == 0) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "sample_id"
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  # Initialize sample_info_mutated
  if (length(meta(object)) != 0) {
    sample_info_mutated <- meta(object, .fmt = "tbl")
  } else {
    sample_info_mutated <- tibble::tibble(sample_id = sample_id(object))
  }
  
  # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
  #----------------------------------------------------------------------------#
  if (length(needed_abundance_keys) > 0) {
    # Create the grid of sample_id and taxa_id combinations
    long_abun <- tidyr::expand_grid(
      sample_id = sample_id(object),
      taxa_id = taxa_id(object)
    )
    
    # Left join all abundance variables needed
    for (abundance_key in needed_abundance_keys) {
      long_abun <- long_abun %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id,
                                names_to = "taxa_id",
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for (i in seq_along(expressions)) {
    
    expr_vars <- all.vars(expressions[[i]])
    
    if (any(expr_vars %in% c("abun", "rela", "norm"))) {
      
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      sample_info_mutated <- long_abun %>%
        dplyr::left_join(sample_info_mutated, by = "sample_id") %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("taxa_id", needed_abundance_keys))) %>%
        dplyr::distinct() %>%
        dplyr::arrange(match(sample_id, sample_id(object)))
      
    } else {
      
      # EXPRESSION WITHOUT ABBUNDANCES
      #----------------------------------------------#
      sample_info_mutated <- sample_info_mutated %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(match(sample_id, sample_id(object)))
      
    }
  } # End of loop over expressions
  
  meta(object) <- tibble::column_to_rownames(sample_info_mutated, "sample_id")
  return(object)
  
})

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
setMethod("mutate_meta", "mgnetList", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (any(nsample(object) == 0)) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","sample_id")
  }
  
  if (!is.character(.by) || "taxa_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'taxa_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#

  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_sample_mutated_merged <- object %>%
    purrr::map(\(x) {
      if(length(meta(x)) != 0){
        meta(x, .fmt = "tbl")
      } else {
        tibble::tibble(sample_id = sample_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
  #----------------------------------------------------------------------------#
  if(length(needed_abundance_keys)!=0){
    
    long_abun_merged <- purrr::map(object, \(mgnet_obj){
      
      long_abun <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                      taxa_id = taxa_id(mgnet_obj))
      
      for(abundance_key in needed_abundance_keys){
        long_abun <- long_abun %>%
          dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                             tibble::as_tibble(rownames = "sample_id") %>%
                             tidyr::pivot_longer(-sample_id, 
                                                 names_to = "taxa_id", 
                                                 values_to = abundance_key),
                           by = c("sample_id", "taxa_id"))
      }
      return(long_abun)
    }) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      info_sample_mutated_merged <- long_abun_merged %>%
        dplyr::left_join(info_sample_mutated_merged, by = c("mgnet","sample_id")) %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      info_sample_mutated_merged <- info_sample_mutated_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  } # loop expressions
  
  info_sample_mutated_splitted <- info_sample_mutated_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(sample_id, sample_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-"mgnet") %>%
        tibble::column_to_rownames("sample_id")
    })
  
  info_sample_mutated_splitted <- info_sample_mutated_splitted[names(object)]
  meta(object) <- info_sample_mutated_splitted
  return(object)
})


#' Modify and Augment `mgnet` Objects by Transforming the `taxa` Slot
#'
#' This function dynamically manipulates the `taxa` slot within `mgnet` or `mgnetList` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'        The function targets the `taxa` slot, which contains metadata for each taxon.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the taxa slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .by Optional; a character vector specifying the columns to group data by before transformations.
#'        Defaults to 'taxa_id' for `mgnet` objects and c('mgnet', 'taxa_id') for `mgnetList` objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by `.by`. Usage of 'sample_id' as a grouping variable is strictly prohibited to maintain
#'        data consistency and avoid misinterpretation.
#'
#'          ### Keywords in `mgnet` and `mgnetList`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **taxa_id**: An essential identifier used to uniquely reference individual taxon within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnetList` object with updated `taxa` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_taxa,mgnet-method mutate_taxa,mgnetList-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate arrange
#' @importFrom tidyr expand_grid any_of
#' @importFrom tidyselect any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_taxa", function(object, ..., .by) {standardGeneric("mutate_taxa")})

setMethod("mutate_taxa", "mgnet", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are taxa to process
  if (ntaxa(object) == 0) {
    stop("Error: No taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- "taxa_id"
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  # Initialize sample_info_mutated
  if (length(taxa(object)) != 0) {
    taxa_info_mutated <- taxa(object, .fmt = "tbl")
  } else {
    taxa_info_mutated <- tibble::tibble(taxa_id = taxa_id(object))
  }
  
  # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
  #----------------------------------------------------------------------------#
  if (length(needed_abundance_keys) > 0) {
    # Create the grid of sample_id and taxa_id combinations
    long_abun <- tidyr::expand_grid(
      sample_id = sample_id(object),
      taxa_id = taxa_id(object)
    )
    
    # Left join all abundance variables needed
    for (abundance_key in needed_abundance_keys) {
      long_abun <- long_abun %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id,
                                names_to = "taxa_id",
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )
    }
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for (i in seq_along(expressions)) {
    
    expr_vars <- all.vars(expressions[[i]])
    
    if (any(expr_vars %in% c("abun", "rela", "norm"))) {
      
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      taxa_info_mutated <- long_abun %>%
        dplyr::left_join(taxa_info_mutated, by = "taxa_id") %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-any_of(c("sample_id", needed_abundance_keys))) %>%
        dplyr::distinct() %>%
        dplyr::arrange(match(taxa_id, taxa_id(object))) %>%
        dplyr::select(-tidyselect::any_of("comm_id"))
      
    } else {
      
      # EXPRESSION WITHOUT ABBUNDANCES
      #----------------------------------------------#
      taxa_info_mutated <- taxa_info_mutated %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(match(taxa_id, taxa_id(object))) %>%
        dplyr::select(-tidyselect::any_of("comm_id"))
      
    }
  } # End of loop over expressions
  
  taxa(object) <- tibble::column_to_rownames(taxa_info_mutated, "taxa_id")
  return(object)
  
})

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
setMethod("mutate_taxa", "mgnetList", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are taxa to process
  if (any(ntaxa(object) == 0)) {
    stop("Error: No taxa available in the 'mgnet' object.")
  }
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
  
  # Check the variables needed
  #lapply(object, \(x){
  #  validate_required_variables(x, keys_required, "taxa")})
  
  
  # Validate the .by argument
  if (missing(.by)) {
    .by <- c("mgnet","taxa_id")
  }
  
  if (!is.character(.by) || "sample_id" %in% .by) {
    stop("Error: '.by' must be a character vector and cannot include 'sample_id'.")
  }
  
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#  
  info_taxa_mutated_merged <- object %>%
    purrr::map(\(x) {
      if(length(taxa(x)) != 0){
        taxa(x, .fmt = "tbl")
      } else {
        tibble::tibble(taxa_id = taxa_id(x))
      }
    }) %>%
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
    purrr::list_rbind()
  
  # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
  #----------------------------------------------------------------------------#
  if(length(needed_abundance_keys)!=0){
    
    long_abun_merged <- purrr::map(object, \(mgnet_obj){
      
      long_abun <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                      taxa_id = taxa_id(mgnet_obj))
      
      for(abundance_key in needed_abundance_keys){
        long_abun <- long_abun %>%
          dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                             tibble::as_tibble(rownames = "sample_id") %>%
                             tidyr::pivot_longer(-sample_id, 
                                                 names_to = "taxa_id", 
                                                 values_to = abundance_key),
                           by = c("sample_id", "taxa_id"))
      }
      return(long_abun)
    }) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
  }
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  for(i in seq_along(expressions)){
    
    expr_vars <- all.vars(expressions[[i]])
    
    if(any(expr_vars %in% c("abun","rela","norm"))){
      
      # EXPRESSION WITH ABUNDANCES
      #----------------------------------------------#
      info_taxa_mutated_merged <- long_abun_merged %>%
        dplyr::left_join(info_taxa_mutated_merged, by = c("mgnet","taxa_id")) %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyr::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() 
      
      
    } else {
      
      # EXPRESSION WITHOUT ABUNDANCES
      #----------------------------------------------#
      info_taxa_mutated_merged <- info_taxa_mutated_merged %>%
        dplyr::group_by(!!!rlang::syms(.by)) %>%
        dplyr::mutate(!!!rlang::eval_tidy(expressions[i])) %>%
        dplyr::ungroup() 
      
    }
  } # loop expressions
  
  info_taxa_mutated_splitted <- info_taxa_mutated_merged %>%
    base::split(.[, "mgnet"]) %>%
    purrr::imap(\(x,y){
      dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
    }) %>%
    purrr::map(\(x){
      x %>% dplyr::select(-any_of(c("mgnet", "comm_id"))) %>%
        tibble::column_to_rownames("taxa_id")
    })
  
  info_taxa_mutated_splitted <- info_taxa_mutated_splitted[names(object)]
  
  taxa(object) <- info_taxa_mutated_splitted
  return(object)
})