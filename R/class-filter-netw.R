#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#' Filter an `mgnet` Object Using Network/Community Slots
#'
#' Filters taxa in an `mgnet` object by evaluating user-defined logical
#' expressions on the network (`netw`) and community (`comm`) slots.
#' Expressions are evaluated in a data mask where the symbols
#' `netw` and `comm` are bound to the current object's slots—or, if the object
#' is grouped via `group_mgnet()`, to the corresponding slots of each subgroup-
#' specific subset.
#'
#' @param object An `mgnet` or `mgnets` object.
#'   This function targets the `netw` (network graph) and `comm` (community/cluster)
#'   slots to compute logical filters, and subsets taxa accordingly.
#' @param ... One or more expressions to be evaluated. Within each expression, the
#'   symbols `netw` and `comm` refer to the network and community slots respectively.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'   from the object at the end.
#' @param .deselect Logical, default `FALSE`. If `TRUE`, remove the link selection
#'   from the object at the end.
#'
#' @details
#' - **Grouping semantics:** This function respects the grouping state set with
#'   `group_mgnet()`. The taxa-level grouping variables are discovered through
#'   `taxa_vars(object)`. If no taxa grouping variables are present, expressions
#'   are evaluated once globally; otherwise they are evaluated *per subgroup*
#'   on a subsetted object, and results are stitched back in the original order.
#' - **No string rewriting:** expressions are evaluated using a `data` mask so the
#'   symbols `netw` and `comm` resolve directly to the intended slots, avoiding
#'   fragile text substitution.
#' - **Selected links:** if `select_link()` has been used earlier, only those edges
#'   are present during the computation (via a temporary subgraph), then the original
#'   network is restored.
#' - **Result type:** each expression must return a logical vector. It may return
#'   either a scalar or a vector of length equal to the number of taxa in scope
#'   (whole object or subgroup).
#'
#' @return The filtered `mgnet` or `mgnets` object.
#'
#' @export
#' @importFrom dplyr mutate select pull filter
#' @importFrom tidyr unite
#' @importFrom rlang enquos get_expr expr_text eval_tidy current_env
#' @importFrom cli cli_abort
#' @importFrom tidyselect all_of
#' @name filter_netw
#' @aliases filter_netw,mgnet-method filter_netw,mgnets-method
setGeneric("filter_netw", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  standardGeneric("filter_netw")
})

setMethod("filter_netw", "mgnet", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  
  # 0) Basic availability checks ------------------------------------------------
  if (miss_slot(object, "netw") && miss_slot(object, "comm")) {
    cli::cli_abort("No network and communities available in the {.cls mgnet} object.")
  }
  
  # 1) Capture user expressions -------------------------------------------------
  quosures <- rlang::enquos(...)
  if (length(quosures) == 0L) {
    cli::cli_abort("{.fn filter_netw} requires at least one expression in `...`.")
  }
  
  # 2) Read taxa grouping from group_mgnet() state ------------------------------
  taxa_groups <- setdiff(get_group_mgnet(object), meta_vars(object))
  if (length(taxa_groups) == 0) taxa_groups <- NULL
  
  # 3) Determine readable names for diagnostics --------------------------------
  expressions_text <- purrr::map_chr(quosures, ~ rlang::expr_text(rlang::get_expr(.x)))
  quosures_names   <- rlang::names2(quosures)
  quosures_names   <- ifelse(quosures_names == "", expressions_text, quosures_names)
  
  # 4) Evaluate expressions -----------------------------------------------------
  filtered_lists <- vector("list", length(quosures))
  
  if (is.null(taxa_groups)) {
    #--------------------------------------------------------------------------#
    #                              UNGROUPED case                              #
    #--------------------------------------------------------------------------#
    mask <- list(
      netw = netw(object, selected = TRUE),
      comm = comm(object)
    )
    
    for (i in seq_along(quosures)) {
      val <- rlang::eval_tidy(quosures[[i]], data = mask, env = rlang::current_env())
      
      if (!(is.logical(val) && is.null(dim(val)))) {
        cli::cli_abort(c(
          "x" = "Expression '{quosures_names[i]}' must return a {.cls logical} vector.",
          "i" = "You supplied a {.cls {class(val)[1]}}{if (!is.null(dim(val))) ' with dimensions' else ''}."
        ))
      }
      
      if (!(length(val) == 1L || length(val) == ntaxa(object))) {
        cli::cli_abort(c(
          "x" = "Expression '{quosures_names[i]}' must return length 1 or {ntaxa(object)}.",
          "i" = "It returned length {length(val)}."
        ))
      }
      
      filtered_lists[[i]] <- if (length(val) == 1L) rep(val, ntaxa(object)) else val
    }
    
  } else {
    #--------------------------------------------------------------------------#
    #                               GROUPED case                               #
    #--------------------------------------------------------------------------#
    subgroups <- taxa(object, .fmt = "tbl") |>
      dplyr::select(tidyselect::all_of(taxa_groups)) |>
      tidyr::unite("_internal_") |>
      dplyr::pull("_internal_")
    
    unique_keys <- unique(subgroups)
    
    for (i in seq_along(quosures)) {
      result <- vector("logical", length = ntaxa(object))
      
      for (key in unique_keys) {
        idx_key <- which(subgroups == key)
        object_subset <- object[, idx_key]
        
        mask_subset <- list(
          netw = netw(object_subset, selected = TRUE),
          comm = comm(object_subset)
        )
        
        val_key <- rlang::eval_tidy(
          quosures[[i]],
          data = mask_subset,
          env = rlang::current_env()
        )
        
        if (!(is.logical(val_key) && is.null(dim(val_key)))) {
          cli::cli_abort(c(
            "x" = "Expression '{quosures_names[i]}' must return a {.cls logical} vector.",
            "i" = "You supplied a {.cls {class(val_key)[1]}}{if (!is.null(dim(val_key))) ' with dimensions' else ''}."
          ))
        }
        
        if (length(val_key) == 1L) {
          result[idx_key] <- rep(val_key, length(idx_key))
        } else if (length(val_key) == length(idx_key)) {
          result[idx_key] <- val_key
        } else {
          cli::cli_abort(
            "Expression '{quosures_names[i]}' returned length {length(val_key)} for a group of size {length(idx_key)}."
          )
        }
      }
      
      filtered_lists[[i]] <- result
    }
  }
  
  # 5) Combine filters with AND semantics --------------------------------------
  filtered_lists <- lapply(filtered_lists, function(v) {
    v[is.na(v)] <- FALSE
    v
  })
  
  final_mask <- purrr::reduce(filtered_lists, `&`)
  new_obj <- object[, final_mask]
  
  # 6) Finalize ----------------------------------------------------------------
  if (isTRUE(.ungroup)) new_obj <- ungroup_mgnet(new_obj)
  if (isTRUE(.deselect)) new_obj <- deselect_link(new_obj)
  validObject(new_obj)
  new_obj
})


setMethod("filter_netw", "mgnets", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  
  if (length(object) == 0L) {
    return(object)
  }
  
  # 0) Basic availability checks ------------------------------------------------
  if (miss_slot(object, "netw", "any") && miss_slot(object, "comm", "any")) {
    cli::cli_abort("No network and communities available in at least one {.cls mgnet} object.")
  }
  
  # 1) Capture once; the single-object method will do the heavy lifting ---------
  quosures <- rlang::enquos(...)
  if (length(quosures) == 0L) {
    cli::cli_abort("{.fn filter_netw} requires at least one expression in `...`.")
  }
  
  # 2) If the collection is grouped, pass the taxa-level groups down to each mgnet
  if (is_mgnet_grouped(object)) {
    out <- purrr::imap(object, function(x, nm) {
      grp <- setdiff(get_group_mgnet(object), meta_vars(object))
      if (length(grp) == 0L) grp <- NULL
      if (!is.null(grp)) {
        x <- group_mgnet(x, !!!rlang::syms(grp))
      }
      filter_netw(x, !!!quosures)
    })
  } else {
    out <- purrr::map(object, ~ filter_netw(.x, !!!quosures))
  }
  
  object@mgnets <- out
  
  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  if (isTRUE(.deselect)) object <- deselect_link(object)
  validObject(object)
  object
})