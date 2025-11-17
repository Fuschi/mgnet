#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#' Modify and Augment `mgnet` Objects by Transforming the `meta` Slot
#'
#' This function dynamically manipulates the `sample` slot within `mgnet` or `mgnets` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnets` object.
#'        The function targets the `sample` slot, which contains metadata for each sample.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the sample slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'        from the object at the end (i.e., drop the grouping attribute).
#'
#' @details The function is designed to integrate seamlessly with the `tidyverse`, allowing users
#'          to utilize familiar and potent data manipulation verbs such as `mutate`, `filter`.
#'          It supports using any `tidyverse`-compatible expressions, including conditional operations,
#'          summarizations, and complex transformations involving both abundance and metadata fields.
#'          This flexibility makes it particularly useful for ecological and biological data analysis,
#'          where combining different data types and conditions is common.
#'
#'          ### Keywords in `mgnet` and `mgnets`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **sample_id**: An essential identifier used to uniquely reference individual samples within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnets` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnets` object with updated `sample` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_meta,mgnet-method mutate_meta,mgnets-method
setGeneric("mutate_meta", function(object, ..., .ungroup = FALSE) {standardGeneric("mutate_meta")})

setMethod("mutate_meta", "mgnet", function(object, ..., .ungroup = FALSE) {
  
  # 1) Quick check for samples
  if (miss_sample(object)) {
    cli::cli_abort(c("x" = "No samples available in the {.cls mgnet} object."))}
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  check_reserved_keywords(exprs)
  #check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for  abundance
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    meta_tbl   <- meta(object, "tbl")
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    taxa_evars <- evars[evars %in% taxa_vars(object)]
    
    if (length(taxa_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(taxa_evars)} taxa variable{?s} ({.val {paste(taxa_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_meta, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "mgnet_mutate_meta_taxa_without_abundance"
      )
    }
  
    # Determine grouping columns
    meta_groups <- setdiff(get_group_mgnet(object), taxa_vars(object))
    if(length(meta_groups) == 0 && is_abun_expr) meta_groups <- "sample_id"
    meta_groups <- rlang::syms(meta_groups)
    
    # Apply mutate logic
    if (is_abun_expr) {
      # Abundance-related expression 
      meta_tbl <- long_abun %>%
        dplyr::left_join(meta_tbl, by = "sample_id")
      if(length(taxa_evars) != 0){
        meta_tbl <- taxa(object, .fmt = "tbl") %>%
          dplyr::select(tidyselect::all_of(c("taxa_id", taxa_evars))) %>%
          dplyr::right_join(meta_tbl, by = "taxa_id")
      }
      meta_tbl <- meta_tbl %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id", "abun", "rela", "norm", taxa_evars))) %>%
        dplyr::distinct()
    } else {
      meta_tbl <- meta_tbl %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
    
    # --------------------------------------------------------------------------
    # 5) Overwrite meta slot and every loop
    # --------------------------------------------------------------------------
    expr_label <- rlang::as_label(exprs[[i]])
    tryCatch(
      {
        meta(object) <- meta_tbl  # this already validates inside the setter
      },
      error = function(err) {
        # Re-throw with context while preserving the original error as parent
        cli::cli_abort(
          c(
            "",
            "x" = "Updating the {.field meta} slot failed while applying {.code {expr_label}}."
          ),
          class  = "mgnet_mutate_meta_assign_error",
          parent = err
        )})
    
  }
  

  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  object
})

#------------------------------------------------------------------------------#
setMethod("mutate_meta", "mgnets", function(object, ..., .ungroup = FALSE) {
  
  # 1) Quick check for samples
  if (miss_sample(object, .fmt = "any")) {
    cli::cli_abort(c("x" = "No samples available in at least one {.cls mgnet} object."))}
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  check_reserved_keywords(exprs)
  #check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for  abundance
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    
    meta_tbl   <- meta(object, .collapse = TRUE)
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    taxa_evars <- evars[evars %in% taxa_vars(object, .fmt = "unique")]
    
    if (length(taxa_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(taxa_evars)} taxa variable{?s} ({.val {paste(taxa_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_meta, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "mgnet_mutate_meta_taxa_without_abundance"
      )
    }
    
    # Determine grouping columns
    meta_groups <- setdiff(get_group_mgnet(object), taxa_vars(object))
    if(length(meta_groups) == 0 && is_abun_expr) meta_groups <- c("mgnet", "sample_id")
    meta_groups <- rlang::syms(meta_groups)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      meta_tbl <- long_abun %>%
        dplyr::left_join(meta_tbl, by = dplyr::join_by(mgnet,sample_id))
      if(length(taxa_evars) != 0){
        meta_tbl <- taxa(object, .collapse = TRUE) %>%
          dplyr::select(tidyselect::all_of(c("mgnet", "taxa_id", taxa_evars))) %>%
          dplyr::right_join(meta_tbl, by = dplyr::join_by(mgnet, taxa_id))
      }
      meta_tbl <- meta_tbl %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id", "abun", "rela", "norm", taxa_evars))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      meta_tbl <- meta_tbl %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
    
    # --------------------------------------------------------------------------
    # 5) Overwrite meta slot and every loop
    # --------------------------------------------------------------------------
    expr_label <- rlang::as_label(exprs[[i]])
    tryCatch(
      {
        meta(object) <- meta_tbl  # this already validates inside the setter
      },
      error = function(err) {
        # Re-throw with context while preserving the original error as parent
        cli::cli_abort(
          c(
            "",
            "x" = "Updating the {.field meta} slot failed while applying {.code {expr_label}}."
          ),
          class  = "mgnet_mutate_meta_assign_error",
          parent = err
        )})
    
  }

  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  object
})



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Modify and Augment `mgnet` Objects by Transforming the `taxa` Slot
#'
#' This function dynamically manipulates the `taxa` slot within `mgnet` or `mgnets` objects,
#' applying user-defined transformations. It leverages the full suite of `tidyverse` tools, particularly
#' `dplyr`, to enable powerful and flexible data transformations.
#'
#' @param object An `mgnet` or `mgnets` object.
#'        The function targets the `taxa` slot, which contains metadata for each taxon.
#' @param ... Dynamic expressions or functions to be applied to the data.
#'        These expressions can manipulate both abundance data (e.g., 'abun', 'rela', 'norm') and
#'        metadata within the taxa slot. This allows for a comprehensive data transformation
#'        experience that supports all standard and custom `tidyverse` manipulation techniques.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'        from the object at the end (i.e., drop the grouping attribute).
#'        
#'          ### Keywords in `mgnet` and `mgnets`:
#'          - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
#'            directly manipulated or used in conjunction with metadata to perform advanced analyses.
#'          - **taxa_id**: An essential identifier used to uniquely reference individual taxon within an `mgnet` object. 
#'          - **mgnet**: Used exclusively within `mgnets` objects to differentiate between multiple `mgnet` objects 
#'            contained in the list.
#'            
#' @return Returns the `mgnet` or `mgnets` object with updated `taxa` slots reflecting the applied transformations.
#'         All other structures within the object remain unchanged, ensuring that only the targeted metadata is modified.
#'
#' @export
#' @aliases mutate_taxa,mgnet-method mutate_taxa,mgnets-method
#' @importFrom dplyr mutate group_by ungroup distinct relocate arrange
#' @importFrom tidyr expand_grid
#' @importFrom tidyselect any_of
#' @importFrom rlang enquos syms quo_get_expr eval_tidy
#' @importFrom purrr map imap list_rbind
#' @importFrom methods slot
#' @importFrom tibble column_to_rownames tibble add_column
setGeneric("mutate_taxa", function(object, ..., .ungroup = FALSE) {standardGeneric("mutate_taxa")})

setMethod("mutate_taxa", "mgnet", function(object, ..., .ungroup = FALSE) {
  
  # 1) Quick check for samples
  if (miss_taxa(object)) {
    cli::cli_abort(c("x" = "No samples available in the {.cls mgnet} object."))}
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  check_reserved_keywords(exprs)
  #check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for  abundance
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    taxa_tbl   <- taxa(object, "tbl")
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    meta_evars <- evars[evars %in% meta_vars(object)]
    
    if (length(meta_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(meta_evars)} taxa variable{?s} ({.val {paste(meta_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_taxa, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "mgnet_mutate_taxa_taxa_without_abundance"
      )
    }
    
    # Determine grouping columns
    taxa_groups <- setdiff(get_group_mgnet(object), meta_vars(object))
    if(length(taxa_groups) == 0 && is_abun_expr) taxa_groups <- "taxa_id"
    taxa_groups <- rlang::syms(taxa_groups)
    
    # Apply mutate logic
    if (is_abun_expr) {
      # Abundance-related expression 
      taxa_tbl <- long_abun %>%
        dplyr::left_join(taxa_tbl, by = "taxa_id")
      if(length(meta_evars) != 0){
        taxa_tbl <- meta(object, .fmt = "tbl") %>%
          dplyr::select(tidyselect::all_of(c("sample_id", meta_evars))) %>%
          dplyr::right_join(taxa_tbl, by = "sample_id")
      }
      taxa_tbl <- taxa_tbl %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm", meta_evars))) %>%
        dplyr::distinct()
    } else {
      taxa_tbl <- taxa_tbl %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
    
    # --------------------------------------------------------------------------
    # 5) Overwrite meta slot and every loop
    # --------------------------------------------------------------------------
    expr_label <- rlang::as_label(exprs[[i]])
    tryCatch(
      {
        taxa(object) <- taxa_tbl 
      },
      error = function(err) {
        # Re-throw with context while preserving the original error as parent
        cli::cli_abort(
          c(
            "",
            "x" = "Updating the {.field taxa} slot failed while applying {.code {expr_label}}."
          ),
          class  = "mgnet_mutate_taxa_assign_error",
          parent = err
        )})
    
  }
  
  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  object
})

#------------------------------------------------------------------------------#
setMethod("mutate_taxa", "mgnets", function(object, ..., .ungroup = FALSE) {
  
  # 1) Quick check for samples
  if (miss_taxa(object, .fmt = "any")) {
    cli::cli_abort(c("x" = "No taxa available in at least one {.cls mgnet} object."))}
  
  # 2) Capture main expressions
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  check_reserved_keywords(exprs)
  #check_forbidden_expressions(exprs)
  
  # 3) Gather data frames for  abundance
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  
  # --------------------------------------------------------------------------
  # 4) Loop over each captured expression
  # --------------------------------------------------------------------------
  for (i in seq_along(exprs)) {
    
    taxa_tbl   <- taxa(object, .collapse = TRUE)
    evars <- all.vars(exprs[[i]])  # variables used in current expression
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    meta_evars <- evars[evars %in% meta_vars(object, .fmt = "unique")]
    
    if (length(meta_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(meta_evars)} taxa variable{?s} ({.val {paste(meta_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_taxa, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "mgnet_mutate_taxa_taxa_without_abundance"
      )
    }
    
    # Determine grouping columns
    taxa_groups <- setdiff(get_group_mgnet(object), meta_vars(object))
    if(length(taxa_groups) == 0 && is_abun_expr) taxa_groups <- c("mgnet", "taxa_id")
    taxa_groups <- rlang::syms(taxa_groups)
    
    # Apply mutate logic
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Abundance-related expression => join abundance first
      taxa_tbl <- long_abun %>%
        dplyr::left_join(taxa_tbl, by = dplyr::join_by(mgnet,taxa_id))
      if(length(meta_evars) != 0){
        taxa_tbl <- meta(object, .collapse = TRUE) %>%
          dplyr::select(tidyselect::all_of(c("mgnet", "sample_id", meta_evars))) %>%
          dplyr::right_join(taxa_tbl, by = dplyr::join_by(mgnet, taxa_id))
      }
      taxa_tbl <- taxa_tbl %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm", meta_evars))) %>%
        dplyr::distinct()
    } else {
      # Non-abundance expression => mutate in-place
      taxa_tbl <- taxa_tbl %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::mutate(!!!rlang::eval_tidy(exprs[i])) %>%
        dplyr::ungroup()
    }
    
    # --------------------------------------------------------------------------
    # 5) Overwrite meta slot and every loop
    # --------------------------------------------------------------------------
    expr_label <- rlang::as_label(exprs[[i]])
    tryCatch(
      {
        taxa(object) <- taxa_tbl  # this already validates inside the setter
      },
      error = function(err) {
        # Re-throw with context while preserving the original error as parent
        cli::cli_abort(
          c(
            "",
            "x" = "Updating the {.field meta} slot failed while applying {.code {expr_label}}."
          ),
          class  = "mgnet_mutate_taxa_assign_error",
          parent = err
        )})
    
  }
  
  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  object
})