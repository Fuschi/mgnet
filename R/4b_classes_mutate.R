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
  
  # 1) Quick check for samples
  if (miss_sample(object)) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for meta & abundance
  meta_df   <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Handle `.by`: 
  #    - If missing => user didn't specify .by
  #    - If NULL    => no grouping
  #    - If character => group by those columns
  # --------------------------------------------------------------------------
  by_is_missing <- missing(.by)
  if (!by_is_missing && !is.null(.by) && !is.character(.by)) {
    stop("`.by` must be either missing, NULL, or a character vector of column names.")
  }
  
  # A helper function to decide the grouping columns per expression
  get_local_group_cols <- function(evars) {
    
    # 1) If `.by` was provided (not missing)
    if (!by_is_missing) {
      # .by is either NULL or character
      if (is.null(.by) || length(.by) == 0) {
        # .by = NULL or .by = character(0) => explicitly no grouping
        return(NULL)
      } else {
        # .by is a non-empty character vector
        return(rlang::syms(.by))
      }
    }
    
    # 2) .by is missing => rely on existing meta_groups or fallback logic
    if (identical(get_group_meta(object), character(0))) {
      # `is_meta_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(get_group_meta(object))) {
      # Use stored meta group columns
      return(rlang::syms(get_group_meta(object)))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms("sample_id"))
    } else {
      # No grouping otherwise
      return(NULL)
    }
  }
  
  # --------------------------------------------------------------------------
  # 5) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    
    # Determine grouping columns
    local_group_cols <- get_local_group_cols(evars)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      meta_df <- long_abun %>%
        dplyr::left_join(meta_df, by = "sample_id") %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      meta_df <- meta_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
  }
  
  # --------------------------------------------------------------------------
  # 6) Overwrite meta slot and return
  # --------------------------------------------------------------------------
  meta(object) <- meta_df
  object
})


setMethod("mutate_meta", "mgnetList", function(object, ..., .by) {
  
  # 1) Quick check for samples
  if (miss_sample(object, "any")) {
    stop("Error: No samples available in at least one 'mgnet' object.")
  }
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for meta & abundance
  meta_df   <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Handle `.by`: 
  #    - If missing => user didn't specify .by
  #    - If NULL    => no grouping
  #    - If character => group by those columns
  # --------------------------------------------------------------------------
  by_is_missing <- missing(.by)
  if (!by_is_missing && !is.null(.by) && !is.character(.by)) {
    stop("`.by` must be either missing, NULL, or a character vector of column names.")
  }
  
  # A helper function to decide the grouping columns per expression
  get_local_group_cols <- function(evars) {
    
    # 1) If `.by` was provided (not missing)
    if (!by_is_missing) {
      # .by is either NULL or character
      if (is.null(.by) || length(.by) == 0) {
        # .by = NULL or .by = character(0) => explicitly no grouping
        return(NULL)
      } else {
        # .by is a non-empty character vector
        return(rlang::syms(.by))
      }
    }
    
    # 2) .by is missing => rely on existing meta_groups or fallback logic
    if (identical(get_group_meta(object), character(0))) {
      # `get_group_meta(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(get_group_meta(object))) {
      # Use stored meta group columns
      return(rlang::syms(get_group_meta(object)))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms(c("mgnet", "sample_id")))
    } else {
      # No grouping otherwise
      return(rlang::syms("mgnet"))
    }
  }
  
  # --------------------------------------------------------------------------
  # 5) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    
    # Determine grouping columns
    local_group_cols <- get_local_group_cols(evars)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      meta_df <- long_abun %>%
        dplyr::left_join(meta_df, by = c("mgnet","sample_id")) %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      meta_df <- meta_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
  }
  
  # --------------------------------------------------------------------------
  # 6) Overwrite meta slot and return
  # --------------------------------------------------------------------------
  meta(object) <- meta_df
  object
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
setGeneric("mutate_taxa", function(object, ..., .by) {standardGeneric("mutate_taxa")})

setMethod("mutate_taxa", "mgnet", function(object, ..., .by) {
  
  # 1) Quick check for samples
  if (miss_taxa(object)) {
    stop("Error: No taxa available in the 'mgnet' object.")
  }
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for taxa & abundance
  taxa_df   <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Handle `.by`: 
  #    - If missing => user didn't specify .by
  #    - If NULL    => no grouping
  #    - If character => group by those columns
  # --------------------------------------------------------------------------
  by_is_missing <- missing(.by)
  if (!by_is_missing && !is.null(.by) && !is.character(.by)) {
    stop("`.by` must be either missing, NULL, or a character vector of column names.")
  }
  
  # A helper function to decide the grouping columns per expression
  get_local_group_cols <- function(evars) {
    
    # 1) If `.by` was provided (not missing)
    if (!by_is_missing) {
      # .by is either NULL or character
      if (is.null(.by) || length(.by) == 0) {
        # .by = NULL or .by = character(0) => explicitly no grouping
        return(NULL)
      } else {
        # .by is a non-empty character vector
        return(rlang::syms(.by))
      }
    }
    
    # 2) .by is missing => rely on existing taxa_groups or fallback logic
    if (identical(get_group_taxa(object), character(0))) {
      # `is_taxa_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(get_group_taxa(object))) {
      # Use stored taxa group columns
      return(rlang::syms(get_group_taxa(object)))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms("taxa_id"))
    } else {
      # No grouping otherwise
      return(NULL)
    }
  }
  
  # --------------------------------------------------------------------------
  # 5) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    
    # Determine grouping columns
    local_group_cols <- get_local_group_cols(evars)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      taxa_df <- long_abun %>%
        dplyr::left_join(taxa_df, by = "taxa_id") %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      taxa_df <- taxa_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
  }
  
  # --------------------------------------------------------------------------
  # 6) Overwrite taxa_df slot and return
  # --------------------------------------------------------------------------
  taxa(object) <- taxa_df
  object
  
})

setMethod("mutate_taxa", "mgnetList", function(object, ..., .by) {
  
  # 1) Quick check for samples
  if (miss_taxa(object, "any")) {
    stop("Error: No taxa available in at least one 'mgnet' object.")
  }
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for taxa & abundance
  taxa_df   <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Handle `.by`: 
  #    - If missing => user didn't specify .by
  #    - If NULL    => no grouping
  #    - If character => group by those columns
  # --------------------------------------------------------------------------
  by_is_missing <- missing(.by)
  if (!by_is_missing && !is.null(.by) && !is.character(.by)) {
    stop("`.by` must be either missing, NULL, or a character vector of column names.")
  }
  
  # A helper function to decide the grouping columns per expression
  get_local_group_cols <- function(evars) {
    
    # 1) If `.by` was provided (not missing)
    if (!by_is_missing) {
      # .by is either NULL or character
      if (is.null(.by) || length(.by) == 0) {
        # .by = NULL or .by = character(0) => explicitly no grouping
        return(NULL)
      } else {
        # .by is a non-empty character vector
        return(rlang::syms(.by))
      }
    }
    
    # 2) .by is missing => rely on existing taxa_groups or fallback logic
    if (identical(get_group_taxa(object), character(0))) {
      # `is_taxa_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(get_group_taxa(object))) {
      # Use stored taxa group columns
      return(rlang::syms(get_group_taxa(object)))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms(c("mgnet","taxa_id")))
    } else {
      # No grouping otherwise
      return("mgnet")
    }
  }
  
  # --------------------------------------------------------------------------
  # 5) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    
    # Determine grouping columns
    local_group_cols <- get_local_group_cols(evars)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      taxa_df <- long_abun %>%
        dplyr::left_join(taxa_df, by = "taxa_id") %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      taxa_df <- taxa_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
  }
  
  # --------------------------------------------------------------------------
  # 6) Overwrite taxa slot and return
  # --------------------------------------------------------------------------
  taxa(object) <- taxa_df
  object
  
})