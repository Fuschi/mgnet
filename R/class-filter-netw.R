#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#' Filter using Network/Community Slots in an `mgnet` Object
#'
#' Applies user-defined expressions to compute new columns for the `taxa` slot,
#' using the network (`netw`) and community (`comm`) slots of an `mgnet` object.
#' Expressions are evaluated in a data mask where the symbols
#' `netw` and `comm` are bound to the current object's slots—or, if the object
#' is grouped via `group_mgnet()`, to the corresponding slots of each subgroup-
#' specific subset.
#'
#' @param object An `mgnet` object.
#'   This function targets the `netw` (network graph) and `comm` (community/cluster)
#'   slots for computations, and filter the object.
#' @param ... One or more expressions to be evaluated. Within each expression, the
#'   symbols `netw` and `comm` refer to the network and community slots respectively.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'        from the object at the end.
#' @param .deselect Logical, default `FALSE`. If `TRUE`, remove the link selection
#'        from the object at the end.
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
#' - **Result length:** each expression must return a vector of
#'   length equal to the number of taxa in scope (whole object or subgroup).
#'
#' @return The modified `mgnet` object, with new columns appended to `taxa(object)`.
#'
#' @export
#' @name filter_netw
#' @aliases filter_netw,mgnet-method filter_netw,mgnets-method
setGeneric("filter_netw", function(object, ..., .ungroup = FALSE, .deselect = FALSE) standardGeneric("filter_netw"))

setMethod("filter_netw", "mgnet", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {

  # 0) Basic availability checks ------------------------------------------------
  if (miss_slot(object, "netw") && miss_slot(object, "comm")) {
    cli::cli_abort("No network and communities available in the {.cls mgnet} object.")
  }

  # 1) Capture user expressions and validate reserved keywords -----------------
  quosures <- rlang::enquos(...)
  if (length(quosures) == 0L) cli::cli_abort("{.fn mutate_netw} requires at least one expression in `...`.")
  check_reserved_keywords(quosures)

  # 2) Read taxa grouping from group_mgnet() state ------------------------------
  #    We assume taxa_vars(object) returns the taxa-level grouping variables
  taxa_groups <- setdiff(get_group_mgnet(object), meta_vars(object))
  if (length(taxa_groups) == 0) taxa_groups <- NULL

  # 3) Determine names for the new output columns ------------------------------
  expressions_text <- purrr::map_chr(quosures, ~ rlang::expr_text(rlang::get_expr(.x)))
  quosures_names  <- rlang::names2(quosures)
  quosures_names  <- ifelse(quosures_names == "", expressions_text, quosures_names)

  # 4) Evaluate expressions (no string replacement; use a data mask) -----------
  filtered_lists <- vector("list", length(quosures))
  if (is.null(taxa_groups)) {
    #--------------------------------------------------------------------------#
    #                             UNGROUPED case                               #
    #--------------------------------------------------------------------------#
    mask <- list(
      netw = netw(object, selected = TRUE),
      comm = comm(object)
    )

    for (i in seq_along(quosures)) {

      val <- rlang::eval_tidy(quosures[[i]], data = mask, env = rlang::current_env())

      if (!(is.logical(val) && is.null(dim(val)))) {
        cli::cli_abort(c(
          "x" = "Expression '{exp_nm}' must return a {.cls logical} vector.",
          "i" = "You supplied a {.cls {class(val)[1]}}{if (!is.null(dim(val))) ' with dimensions' else ''}."))}

      if (length(val) != ntaxa(object)) {
        cli::cli_abort(
          "Expression '{exp_nm}' must return a logical vector of length {n}, not {length(val)}.")}

      filtered_lists[[i]] <- val
    }

  } else {
    #--------------------------------------------------------------------------#
    #                             UNGROUPED case                               #
    #--------------------------------------------------------------------------#
    # Build subgroup keys from taxa metadata
    subgroups <- taxa(object, .fmt = "tbl") |>
      dplyr::select(tidyselect::all_of(taxa_groups)) |>
      tidyr::unite("_internal_") |>
      dplyr::pull("_internal_")

    unique_keys <- unique(subgroups)

    for (i in seq_along(quosures)) {
      # Store per-taxa results; allow per-group scalars or vectors
      result <- vector(length = ntaxa(object))

      for (key in unique_keys) {
        idx_key <- which(subgroups == key)
        object_subset <- object[, idx_key]

        # Per-group mask: bind netw/comm of the subset
        mask_subset <- list(
          netw = netw(object_subset, selected = TRUE),
          comm = comm(object_subset)
        )

        val_key <- rlang::eval_tidy(quosures[[i]], data = mask_subset, env = rlang::current_env())

        # Recycle scalar or validate vector length for this group
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
      
      if (!(is.logical(result) && is.null(dim(result)))) {
        cli::cli_abort(c(
          "x" = "Expression '{exp_nm}' must return a {.cls logical} vector.",
          "i" = "You supplied a {.cls {class(result)[1]}}{if (!is.null(dim(result))) ' with dimensions' else ''}."))}
      
      if (length(result) != ntaxa(object)) {
        cli::cli_abort(
          "Expression '{exp_nm}' must return a logical vector of length {ntaxa(object)}, not {length(result)}.")}
      

      filtered_lists[[i]] <- result
    }
  }

  # 5) Intersect sample sets from all expressions (AND logic)
  filtered_lists <- lapply(filtered_lists, function(v) { v[is.na(v)] <- FALSE; v }) # Replace NA to FALSE as dplyr style
  final_mask <- purrr::reduce(filtered_lists, `&`)   # AND for each element
  new_obj <- object[, final_mask]

  # 5) Validate and return ------------------------------------------------------
  if(isTRUE(.ungroup)) new_obj <- ungroup_mgnet(new_obj)
  if(isTRUE(.deselect)) new_obj <- deselect_link(new_obj)
  new_obj
})


#' @rdname filter_netw
#' @export
setMethod("filter_netw", "mgnets", function(object, ..., .ungroup = FALSE, .deselect = FALSE) {
  # Early return if there are no elements
  if (length(object@mgnets) == 0L) return(object)
  
  # Compute mgnets-level taxa groups once and sanitize them
  grp <- setdiff(get_group_mgnet(object), meta_vars(object))
  if (length(grp) == 0L) grp <- NULL
  
  # Apply filter_netw to each mgnet element.
  object@mgnets <- purrr::imap(object@mgnets, function(x, nm) {
    
    # Apply groups from the mgnets object (if any)
    grp <- setdiff(get_group_mgnet(object), meta_vars(object))
    if (length(grp) == 0L) grp <- NULL
    if (!is.null(grp)) {
      x <- tryCatch(
        rlang::inject(group_mgnet(x, !!!rlang::syms(grp))),
        error = function(e) cli::cli_abort(c(
          "Failed to apply groups to element {label}.",
          "x" = e$message,
          "i" = "Ensure {.fn group_mgnet} for {.cls mgnet} accepts tidy-select columns."
        ))
      )
    }
    
    # Run filter
    tryCatch(
      filter_netw(x, ..., .ungroup = FALSE, .deselect = FALSE),
      error = function(e) cli::cli_abort(c(
        "Failed while filtering element {label} of the {.cls mgnets} object.",
        "x" = e$message,
        "i" = "Each expression in `...` must evaluate against {.code netw}/{.code comm} and return a logical vector."
      ))
    )
  })
  
  # Apply ungroup/deselect at the mgnets level (once, at the very end)
  if (isTRUE(.ungroup))  object <- ungroup_mgnet(object)
  if (isTRUE(.deselect)) object <- deselect_link(object)
  
  object
})
