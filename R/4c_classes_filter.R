#' Filter `mgnet` or `mgnetList` Objects by Sample Metadata and Abundance Fields
#'
#' Filters the samples in an `mgnet` or `mgnetList` object based on user-specified 
#' conditions provided through `...`. Each condition can reference columns in the 
#' sample metadata (`meta`) or abundance-related fields (`abun`, `rela`, `norm`). 
#' When multiple conditions are provided, they are combined using a logical AND, 
#' meaning a sample must satisfy **all** conditions to be retained.
#'
#' @param object An `mgnet` or `mgnetList` object. For `mgnetList`, each contained 
#'   `mgnet` can be filtered according to the combined conditions, with a column 
#'   `mgnet` distinguishing them in any joined data frames.
#' @param ... Unquoted expressions for filtering. These expressions can reference 
#'   columns in the sample metadata (`meta`) or any of the abundance-related fields 
#'   (`abun`, `rela`, `norm`). Each expression is applied sequentially, and only 
#'   samples that satisfy **all** expressions are retained.
#' @param .by A character vector of grouping columns. The default differs by object 
#'   type:
#'   \itemize{
#'     \item For `mgnet`, the default is `.by = "sample_id"`.
#'     \item For `mgnetList`, the default is `.by = c("mgnet", "sample_id")`.
#'   }
#'   If `.by` is `NULL`, no grouping is applied, and conditions are evaluated 
#'   across all samples globally. If `.by` is missing, any existing meta grouping 
#'   (from `group_meta()`) is used instead.
#'
#' @details
#' This function integrates closely with **tidyverse** conventions. Each expression 
#' in `...` is captured and evaluated using tidy evaluation, allowing you to write 
#' conditions similarly to `dplyr::filter()`. When an expression references 
#' abundance variables (`abun`, `rela`, or `norm`), the function automatically 
#' joins the sample metadata with the relevant abundance data (`long_abun`) before 
#' filtering.
#'
#' **Grouping Logic**:
#' 1. **If you provide an explicit `.by`** (a non-empty character vector), those 
#'    columns define the grouping.
#' 2. **If `.by` is `NULL`**, no grouping is applied.
#' 3. **If `.by` is missing**:
#'    - **If the object has existing meta grouping** (set via `group_meta()`), that 
#'      grouping is used.
#'    - **Else**, if any filter expression references abundance variables, 
#'      grouping defaults to:
#'      \itemize{
#'        \item `c("mgnet", "sample_id")` for `mgnetList` objects.
#'        \item `"sample_id"` for `mgnet` objects.
#'      }
#'    - **Otherwise**, no grouping is applied.
#'
#' **Multiple Conditions**:
#' Each condition in `...` is applied sequentially, narrowing down the set of 
#' remaining samples. After processing all conditions, only those samples that 
#' satisfied **all** conditions remain.
#'
#' **Behavior with `mgnetList`**:
#' - The `mgnet` column in the joined data frames identifies which `mgnet` each sample 
#'   belongs to.
#' - Filtering expressions can apply to metadata or abundance columns, referencing 
#'   each `mgnet`'s data.
#' - Samples are filtered **independently** for each `mgnet` based on the conditions.
#'
#' **Return Value**:
#' - For an `mgnet` object, an updated `mgnet` containing only the samples that 
#'   passed all filters.
#' - For a `mgnetList`, an updated list of `mgnet` objects. Each `mgnet` is 
#'   subset to the remaining samples.
#'
#' This approach provides a flexible, `dplyr`-like interface for complex filtering 
#' of both sample metadata and abundance data, respecting the specified grouping context.
#'
#' @return An `mgnet` or `mgnetList` object containing only the samples (and their 
#'   corresponding abundance data, if any) that satisfy all provided filter conditions. 
#'   The rest are removed. If no samples match, an object with zero samples is returned 
#'   (but the same structure otherwise).
#'
#' @export
#' @aliases filter_meta,mgnet-method filter_meta,mgnetList-method
setGeneric("filter_meta", function(object, ..., .by) {
  standardGeneric("filter_meta")
})

#' @rdname filter_meta
#' @export
setMethod("filter_meta", "mgnet", function(object, ..., .by) {
  
  # 1) Check for samples
  if (miss_sample(object)) {
    stop("Error: No samples available in the 'mgnet' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  meta_df   <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  if(length(long_abun) && length(meta_df)){
    long_abun <- dplyr::left_join(long_abun, meta_df, by = "sample_id") 
  }
  
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
    stored_groups <- get_group_meta(object)
    if (identical(stored_groups, character(0))) {
      # `is_meta_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(stored_groups)) {
      # Use stored meta group columns
      return(rlang::syms(stored_groups))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms("sample_id"))
    } else {
      # No grouping otherwise
      return(NULL)
    }
  }
  
  # 5) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    
    # Determine grouping
    local_group_cols <- get_local_group_cols(evars)
    
    # Decide if expression references abun/rela/norm
    if (any(evars %in% c("abun","rela","norm"))) {
      # Join abundance
      filtered_lists[[i]] <- long_abun %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id","abun","rela","norm"))) %>%
        dplyr::distinct() %>%
        dplyr::pull("sample_id")
    } else {
      filtered_lists[[i]] <- meta_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::pull("sample_id")
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_samples <- purrr::reduce(filtered_lists, intersect)
  
  # Reorder final_samples to match the original object’s sample order
  sample_ids_in_object <- sample_id(object)
  final_samples <- final_samples[order(match(final_samples, sample_ids_in_object))]
  
  # 7) Subset the object by the filtered samples and return
  return(object[final_samples, ])
})

#' @rdname filter_meta
#' @export
setMethod("filter_meta", "mgnetList", function(object, ..., .by) {
  
  # 1) Check for samples
  if (miss_sample(object, "any")) {
    stop("Error: No samples available in at least one 'mgnet' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  meta_df   <- gather_meta(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  if(length(long_abun) && length(get_abundance_keys(exprs))){
    long_abun <- dplyr::left_join(long_abun, meta_df, by = c("mgnet", "sample_id")) 
  }
  
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
    stored_groups <- get_group_meta(object)
    if (identical(stored_groups, character(0))) {
      # `is_meta_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(stored_groups)) {
      # Use stored meta group columns
      return(rlang::syms(stored_groups))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to sample_id if abundance variables are used
      return(rlang::syms(c("mgnet", "sample_id")))
    } else {
      # No grouping otherwise
      return(rlang::syms("mgnet"))
    }
  }
  
  # 5) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    
    # Determine grouping
    local_group_cols <- get_local_group_cols(evars)
    # print(get_local_group_cols(evars))  # Uncomment if debugging is needed
    
    # Decide if expression references abun/rela/norm
    if (any(evars %in% c("abun","rela","norm"))) {
      # Join abundance data with metadata
      filtered_lists[[i]] <- long_abun %>%
        dplyr::left_join(meta_df, by = c("mgnet","sample_id")) %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id","abun","rela","norm"))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyselect::all_of(c("mgnet", "sample_id")))
    } else {
      # Filter metadata alone
      filtered_lists[[i]] <- meta_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(tidyselect::all_of(c("mgnet", "sample_id")))
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_samples <- purrr::reduce(filtered_lists, dplyr::semi_join, by = c("mgnet", "sample_id")) 
  
  # 7) Subset each mgnet in the mgnetList by the final filtered samples
  for(i in names(object)){
    # Extract the sample_ids for this mgnet that passed all filters
    filtered_i <- dplyr::filter(final_samples, mgnet == i) %>% dplyr::pull("sample_id")
    
    # Reorder survived_ids to match the original sample order
    object_samples <- sample_id(object[[i]])
    matched_indices <- match(filtered_i, object_samples)
    matched_indices <- matched_indices[!is.na(matched_indices)]
    filtered_i <- object_samples[matched_indices]
    object[[i]] <- object[[i]][filtered_i, ]
  }
  
  # 8) Validate and return the updated mgnetList
  validObject(object)
  return(object)
})






#' Filter `mgnet` or `mgnetList` Objects by Taxa Metadata and Abundance Fields
#'
#' Filters the taxa in an `mgnet` or `mgnetList` object based on user-specified 
#' conditions provided through `...`. Each condition can reference columns in the 
#' taxa metadata (`taxa`) or abundance-related fields (`abun`, `rela`, `norm`). 
#' When multiple conditions are provided, they are combined using a logical AND, 
#' meaning a taxon must satisfy **all** conditions to be retained.
#'
#' @param object An `mgnet` or `mgnetList` object. For `mgnetList`, each contained 
#'   `mgnet` can be filtered according to the combined conditions, with a column 
#'   `mgnet` distinguishing them in any joined data frames.
#' @param ... Unquoted expressions for filtering. These expressions can reference 
#'   columns in the taxa metadata (`meta`) or any of the abundance-related fields 
#'   (`abun`, `rela`, `norm`). Each expression is applied sequentially, and only 
#'   taxa that satisfy **all** expressions are retained.
#' @param .by A character vector of grouping columns. The default differs by object 
#'   type:
#'   \itemize{
#'     \item For `mgnet`, the default is `.by = "taxa_id"`.
#'     \item For `mgnetList`, the default is `.by = c("mgnet", "taxa_id")`.
#'   }
#'   If `.by` is `NULL`, no grouping is applied, and conditions are evaluated 
#'   across all taxa globally. If `.by` is missing, any existing taxa grouping 
#'   (from `group_taxa()`) is used instead.
#'
#' @details
#' This function integrates closely with **tidyverse** conventions. Each expression 
#' in `...` is captured and evaluated using tidy evaluation, allowing you to write 
#' conditions similarly to `dplyr::filter()`. When an expression references 
#' abundance variables (`abun`, `rela`, or `norm`), the function automatically 
#' joins the taxa metadata with the relevant abundance data (`long_abun`) before 
#' filtering.
#'
#' **Grouping Logic**:
#' 1. **If you provide an explicit `.by`** (a non-empty character vector), those 
#'    columns define the grouping.
#' 2. **If `.by` is `NULL`**, no grouping is applied.
#' 3. **If `.by` is missing**:
#'    - **If the object has existing taxa grouping** (set via `group_taxa()`), that 
#'      grouping is used.
#'    - **Else**, if any filter expression references abundance variables, 
#'      grouping defaults to:
#'      \itemize{
#'        \item `c("mgnet", "taxa_id")` for `mgnetList` objects.
#'        \item `"taxa_id"` for `mgnet` objects.
#'      }
#'    - **Otherwise**, no grouping is applied.
#'
#' **Multiple Conditions**:
#' Each condition in `...` is applied sequentially, narrowing down the set of 
#' remaining taxa. After processing all conditions, only those taxa that 
#' satisfied **all** conditions remain.
#'
#' **Behavior with `mgnetList`**:
#' - The `mgnet` column in the joined data frames identifies which `mgnet` each taxon 
#'   belongs to.
#' - Filtering expressions can apply to metadata or abundance columns, referencing 
#'   each `mgnet`'s data.
#' - Taxa are filtered **independently** for each `mgnet` based on the conditions.
#'
#' **Return Value**:
#' - For an `mgnet` object, an updated `mgnet` containing only the taxa that 
#'   passed all filters.
#' - For a `mgnetList`, an updated list of `mgnet` objects. Each `mgnet` is 
#'   subset to the remaining taxa.
#'
#' This approach provides a flexible, `dplyr`-like interface for complex filtering 
#' of both taxa metadata and abundance data, respecting the specified grouping context.
#'
#' @return An `mgnet` or `mgnetList` object containing only the taxa (and their 
#'   corresponding abundance data, if any) that satisfy all provided filter conditions. 
#'   The rest are removed. If no taxa match, an object with zero taxa is returned 
#'   (but the same structure otherwise).
#'
#' @export
#' @aliases filter_taxa,mgnet-method filter_taxa,mgnetList-method
setGeneric("filter_taxa", function(object, ..., .by) {
  standardGeneric("filter_taxa")
})

#' @rdname filter_taxa
#' @export
setMethod("filter_taxa", "mgnet", function(object, ..., .by) {
  
  # 1) Check for taxa
  if (miss_taxa(object)) {
    stop("Error: No taxa available in the 'mgnet' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  taxa_df   <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # Optional: Join abundance with taxa metadata if both are present
  if (length(long_abun) && length(get_abundance_keys(exprs))) {
    long_abun <- dplyr::left_join(long_abun, taxa_df, by = "taxa_id") 
  }
  
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
    stored_groups <- get_group_taxa(object)
    if (identical(stored_groups, character(0))) {
      # `taxa_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(stored_groups)) {
      # Use stored taxa group columns
      return(rlang::syms(stored_groups))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to taxa_id if abundance variables are used
      return(rlang::syms("taxa_id"))
    } else {
      # No grouping otherwise
      return(NULL)
    }
  }
  
  # 5) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    
    # Determine grouping
    local_group_cols <- get_local_group_cols(evars)
    
    # Decide if expression references abun/rela/norm
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Filter using abundance data
      filtered_lists[[i]] <- long_abun %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        # Remove abundance columns to avoid duplication
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm"))) %>%
        dplyr::distinct() %>%
        dplyr::pull("taxa_id")
    } else {
      # Filter using taxa taxadata alone
      filtered_lists[[i]] <- taxa_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::pull("taxa_id")
    }
  }
  
  # 6) Intersect taxa sets from all expressions (AND logic)
  final_taxa <- purrr::reduce(filtered_lists, intersect)
  
  # Reorder final_taxa to match the original object’s taxa order
  taxa_ids_in_object <- taxa_id(object)
  final_taxa <- final_taxa[order(match(final_taxa, taxa_ids_in_object))]
  
  # 7) Subset the object by the filtered taxa and return
  return(object[, final_taxa])
})

#' @rdname filter_taxa
#' @export
setMethod("filter_taxa", "mgnetList", function(object, ..., .by) {
  
  # 1) Check for taxa
  if (miss_taxa(object, "any")) {
    stop("Error: No taxa available in at least one 'mgnet' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  taxa_df   <- gather_taxa(object)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # Optional: Join abundance with taxa taxadata if both are present
  if (length(long_abun) && length(get_abundance_keys(exprs))) {
    long_abun <- dplyr::left_join(long_abun, taxa_df, by = c("mgnet","taxa_id")) 
  }
  
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
    stored_groups <- get_group_taxa(object)
    if (identical(stored_groups, character(0))) {
      # `taxa_grouped(object)` returns character(0) => no grouping
      return(NULL)
    } else if (length(stored_groups)) {
      # Use stored taxa group columns
      return(rlang::syms(stored_groups))
    } else if (any(evars %in% c("abun", "rela", "norm"))) {
      # Fallback to mgnet and taxa_id if abundance variables are used
      return(rlang::syms(c("mgnet", "taxa_id")))
    } else {
      # No grouping otherwise
      return(rlang::syms("mgnet"))
    }
  }
  
  # 5) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    
    # Determine grouping
    local_group_cols <- get_local_group_cols(evars)

    # Decide if expression references abun/rela/norm
    if (any(evars %in% c("abun","rela","norm"))) {
      # Join abundance data with taxadata
      filtered_lists[[i]] <- long_abun %>%
        dplyr::left_join(taxa_df, by = c("mgnet","taxa_id")) %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id","abun","rela","norm"))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyselect::all_of(c("mgnet", "taxa_id")))
    } else {
      # Filter taxadata alone
      filtered_lists[[i]] <- taxa_df %>%
        dplyr::group_by(!!!local_group_cols) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(tidyselect::all_of(c("mgnet", "taxa_id")))
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_taxa <- purrr::reduce(filtered_lists, dplyr::semi_join, by = c("mgnet", "taxa_id")) 
  
  # 7) Subset each mgnet in the mgnetList by the final filtered samples
  for(i in names(object)){
    # Extract the sample_ids for this mgnet that passed all filters
    filtered_i <- dplyr::filter(final_taxa, mgnet == i) %>% dplyr::pull("taxa_id")
    
    # Reorder survived_ids to match the original sample order
    object_taxa <- taxa_id(object[[i]])
    matched_indices <- match(filtered_i, object_taxa)
    matched_indices <- matched_indices[!is.na(matched_indices)]
    filtered_i <- object_taxa[matched_indices]
    object[[i]] <- object[[i]][, filtered_i]
  }
  
  # 8) Validate and return the updated mgnetList
  validObject(object)
  return(object)
})
