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
#'        Defaults to 'sample_id' for mgnet objects and c('mgnet', 'sample_id') for mgnetList objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by .by. If you do not wish to group the data, set .by to NULL.
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
  if (miss_sample(object)) stop("Error: No samples available in the 'mgnet' object.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- "sample_id"
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  meta_mutated <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  meta(object) <- apply_mutate_verb(meta_mutated, long_abun, expressions, .by, "sample")
  return(object)
  
})

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
setMethod("mutate_meta", "mgnetList", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if(miss_sample(object, "any")) stop("Error: No sample available in at least one of the mgnet objects.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- c("mgnet", "sample_id")
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  meta_mutated_merged <- gather_meta(object)
  long_abun_merged <- long_abundance_join(object, get_abundance_keys(expressions)) 
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  meta_mutated_merged <- apply_mutate_verb(meta_mutated_merged, long_abun_merged, expressions, .by, "sample")
  
  meta_mutated_splitted <- split_arrange_merged_meta(meta_mutated_merged, object)
  meta(object) <- meta_mutated_splitted
  return(object)
})


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
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
#'        Defaults to 'taxa_id' for mgnet objects and c('mgnet', 'taxa_id') for mgnetList objects.
#'        Grouping ensures that transformations are contextually applied within each subgroup defined
#'        by .by. If you do not wish to group the data, set .by to NULL.
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
#' @importFrom tidyr expand_grid
#' @importFrom tidyselect any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_taxa", function(object, ..., .by = NULL) {standardGeneric("mutate_taxa")})

setMethod("mutate_taxa", "mgnet", function(object, ..., .by = NULL) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are samples to process
  if (miss_taxa(object)) stop("Error: No taxa available in the 'mgnet' object.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)

  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- "taxa_id"
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  taxa_mutated <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(expressions))  

  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  taxa_mutated <- apply_mutate_verb(taxa_mutated, long_abun, expressions, .by, "taxa") %>%
    dplyr::select(-tidyselect::any_of("comm_id"))
  
  taxa(object) <- tibble::column_to_rownames(taxa_mutated, "taxa_id")
  return(object)
  
})

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
setMethod("mutate_taxa", "mgnetList", function(object, ..., .by) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  # Ensure there are taxa to process
  if(miss_taxa(object, "any")) stop("Error: No taxa available in at least one of the mgnet objects.")
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Initialize groups with default values if it is empty.
  if(missing(.by)) .by <- c("mgnet", "taxa_id")
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  taxa_mutated_merged <- gather_taxa(object)
  long_abun_merged <- long_abundance_join(object, get_abundance_keys(expressions)) 
  
  # LOOP OVER THE EXPRESSIONS
  #----------------------------------------------------------------------------#
  taxa_mutated_merged <- apply_mutate_verb(taxa_mutated_merged, long_abun_merged, expressions, .by, "taxa")

  taxa_mutated_splitted <- split_arrange_merged_taxa(taxa_mutated_merged, object)
  taxa(object) <- taxa_mutated_splitted
  return(object)
})