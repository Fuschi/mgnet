#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

# =============================================================================
# Docs (generics share these @rdname blocks)
# =============================================================================

#' Mutate `meta` or `taxa` slots in `mgnet` / `mgnets`
#'
#' Adds or modifies columns in the sample metadata (`meta`) or taxa metadata
#' (`taxa`) of an `mgnet` object. For `mgnets`, the same operation is applied
#' to each contained `mgnet`.
#'
#' User expressions are evaluated with tidy evaluation (like `dplyr::mutate()`).
#' Expressions may reference columns in the current target table (`meta` or `taxa`).
#' If abundance fields are referenced (`abun`, `rela`, `norm`), the relevant
#' long-format abundance table is automatically joined to provide the necessary
#' context.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param ... One or more tidy-evaluated expressions (as in `dplyr::mutate()`).
#'   Named expressions create/overwrite columns with the given name.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'   from the object at the end.
#'
#' @details
#' ## Grouping semantics
#' These functions respect the grouping state set by `group_mgnet()`.
#' Grouping variables are retrieved via `get_group_mgnet(object)` and used to
#' define the context for `dplyr::group_by()` during evaluation.
#'
#' - If the object is **not grouped** (`get_group_mgnet(object)` is `NULL`),
#'   no grouping is applied.
#' - If an expression **references abundance fields** (`abun`, `rela`, `norm`)
#'   and **no explicit grouping** is set, a minimal grouping is injected:
#'   - for `mgnet`: by `sample_id` (when mutating `meta`) or by `taxa_id`
#'     (when mutating `taxa`);
#'   - for `mgnets`: by `mgnet` + the corresponding id.
#' - If an expression **does not reference abundance fields**, grouping variables
#'   that belong to the *other* slot are ignored.
#'
#' ## Cross-slot variables
#' If an expression references variables from the other slot, the expression
#' must also reference at least one abundance field (`abun`, `rela`, `norm`),
#' because cross-slot information is only available through abundance-aware joins.
#'
#' @return An updated `mgnet` or `mgnets` object.
#' @name mutate_mgnet
#' @aliases mutate_meta mutate_taxa
NULL

#' Filter `meta` or `taxa` in `mgnet` / `mgnets`
#'
#' Filters samples (`filter_meta()`) or taxa (`filter_taxa()`) in an `mgnet` object
#' using one or more tidy-evaluated logical expressions (like `dplyr::filter()`).
#' For `mgnets`, filtering is applied independently to each contained `mgnet`.
#'
#' Expressions may reference columns in the target table (`meta` or `taxa`).
#' If abundance fields are referenced (`abun`, `rela`, `norm`), the relevant
#' long-format abundance table is automatically joined to provide the necessary
#' context.
#'
#' Multiple expressions are combined using logical AND.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param ... One or more tidy-evaluated filtering expressions.
#' @param .ungroup Logical, default `FALSE`. If `TRUE`, remove the mgnet grouping
#'   from the object at the end.
#'
#' @details
#' See `mutate_mgnet` for grouping and cross-slot semantics.
#'
#' @return A filtered `mgnet` or `mgnets` object.
#' @name filter_mgnet
#' @aliases filter_meta filter_taxa
NULL

# =============================================================================
# Context helpers (mgnet vs mgnets + meta vs taxa)
# =============================================================================

.mgnet_ctx <- function(object, target = c("meta", "taxa")) {
  is_list <- methods::is(object, "mgnets")
  target  <- match.arg(target)
  cross   <- if (target == "meta") "taxa" else "meta"
  
  # Target/cross ids
  target_key <- if (target == "meta") "sample_id" else "taxa_id"
  cross_key  <- if (cross  == "meta") "sample_id" else "taxa_id"
  
  # Join keys (include mgnet for lists)
  join_target_key <- if (is_list) c("mgnet", target_key) else target_key
  join_cross_key  <- if (is_list) c("mgnet", cross_key)  else cross_key
  
  # Getters (tibbles)
  if (!is_list) {
    get_target_tbl <- function(x) if (target == "meta") meta(x, .fmt = "tbl", .empty = "id") else taxa(x, .fmt = "tbl", .empty = "id")
    get_cross_tbl  <- function(x) if (cross  == "meta") meta(x, .fmt = "tbl", .empty = "id") else taxa(x, .fmt = "tbl", .empty = "id")
  } else {
    get_target_tbl <- function(x) if (target == "meta") meta(x, .collapse = TRUE, .empty = "id") else taxa(x, .collapse = TRUE, .empty = "id")
    get_cross_tbl  <- function(x) if (cross  == "meta") meta(x, .collapse = TRUE, .empty = "id") else taxa(x, .collapse = TRUE, .empty = "id")
  }
  
  # Variable pools
  get_target_vars <- function(x) if (target == "meta") meta_vars(x, .fmt = "unique") else taxa_vars(x, .fmt = "unique")
  get_cross_vars  <- function(x) if (cross  == "meta") meta_vars(x, .fmt = "unique") else taxa_vars(x, .fmt = "unique")
  
  list(
    target = target, cross = cross, is_list = is_list,
    target_key = target_key, cross_key = cross_key,
    join_target_key = join_target_key, join_cross_key = join_cross_key,
    get_target_tbl = get_target_tbl, get_cross_tbl = get_cross_tbl,
    get_target_vars = get_target_vars, get_cross_vars = get_cross_vars
  )
}

# =============================================================================
# Errors
# =============================================================================

.mgnet_abort_missing_slot <- function(object, ctx) {
  if (!ctx$is_list) {
    if (ctx$target == "meta" && miss_sample(object)) {
      cli::cli_abort(c("x" = "No samples available in the {.cls mgnet} object."),
                     class = c("mgnet_no_samples", "mgnet_error"))
    }
    if (ctx$target == "taxa" && miss_metataxa(object)) {
      cli::cli_abort(c("x" = "No taxa available in the {.cls mgnet} object."),
                     class = c("mgnet_no_taxa", "mgnet_error"))
    }
    return(invisible(NULL))
  }
  
  if (ctx$target == "meta" && miss_sample(object, "any")) {
    cli::cli_abort(c("x" = "No samples available in at least one {.cls mgnet} object."),
                   class = c("mgnet_no_samples_any", "mgnet_error"))
  }
  if (ctx$target == "taxa" && miss_metataxa(object, "any")) {
    cli::cli_abort(c("x" = "No taxa available in at least one {.cls mgnet} object."),
                   class = c("mgnet_no_taxa_any", "mgnet_error"))
  }
  
  invisible(NULL)
}

.mgnet_abort_cross_without_abun <- function(expr, cross_vars, verb, ctx) {
  if (length(cross_vars) == 0L) return(invisible(NULL))
  
  cli::cli_abort(
    c(
      "x" = "You referenced {ctx$cross} variable{?s} ({.val {paste(cross_vars, collapse = ', ')}}) without using abundance.",
      "i" = "When {ctx$cross} variables appear in {.fn {paste0(verb, '_', ctx$target)}}, you must also use at least one of {.field abun}, {.field rela}, {.field norm}.",
      "*" = "Expression: {.code {rlang::as_label(expr)}}"
    ),
    class = c(
      paste0("mgnet_", verb, "_", ctx$target, "_", ctx$cross, "_without_abundance"),
      "mgnet_cross_without_abundance",
      "mgnet_error"
    )
  )
}

.mgnet_abort_assign_failed <- function(err, verb, ctx, expr_label) {
  cli::cli_abort(
    c(
      "x" = "Updating the {.field {ctx$target}} slot failed.",
      "*" = "Call: {.code {paste0(verb, '_', ctx$target)}}",
      "*" = "Expression: {.code {expr_label}}"
    ),
    class  = c(paste0("mgnet_", verb, "_", ctx$target, "_assign_error"), "mgnet_error"),
    parent = err
  )
}

.mgnet_abort_filter_failed <- function(err, ctx, stage, expr_label = NULL) {
  bullets <- c(
    "x" = "Filtering failed in {.fn {paste0('filter_', ctx$target)}}.",
    "*" = "Stage: {.val {stage}}"
  )
  if (!is.null(expr_label)) bullets <- c(bullets, "*" = "Expression: {.code {expr_label}}")
  
  cli::cli_abort(
    bullets,
    class  = c(paste0("mgnet_filter_", ctx$target, "_failed"), "mgnet_error"),
    parent = err
  )
}

# =============================================================================
# Expression analysis + validation
# =============================================================================

.mgnet_info_quo <- function(object, ctx, quo, abun_vars = c("abun", "rela", "norm")) {
  vars_used <- all.vars(rlang::get_expr(quo))
  
  cross_picked <- intersect(vars_used, ctx$get_cross_vars(object))
  abun_picked  <- intersect(vars_used, abun_vars)
  
  list(
    expr = quo,
    cross_vars = cross_picked,
    abun_keys  = abun_picked,
    uses_abun  = length(abun_picked) > 0L
  )
}

.mgnet_info_exprs <- function(object, ctx, exprs, abun_vars = c("abun", "rela", "norm")) {
  if (length(exprs) == 0L) {
    return(list(cross_vars = character(0), abun_keys = character(0), uses_abun = FALSE))
  }
  
  infos <- lapply(exprs, function(q) .mgnet_info_quo(object, ctx, q, abun_vars = abun_vars))
  
  cross_all <- unique(unlist(lapply(infos, `[[`, "cross_vars")))
  abun_all  <- unique(unlist(lapply(infos, `[[`, "abun_keys")))
  
  list(
    cross_vars = cross_all,
    abun_keys  = abun_all,
    uses_abun  = length(abun_all) > 0L
  )
}

.mgnet_validate_exprs <- function(object, ctx, exprs, verb, abun_vars = c("abun", "rela", "norm")) {
  if (length(exprs) == 0L) return(invisible(NULL))
  
  # LHS must not overwrite reserved keywords
  check_reserved_keywords(exprs)
  
  # Cross-vars allowed only if the same quo uses abundance
  for (q in exprs) {
    info <- .mgnet_info_quo(object, ctx, q, abun_vars = abun_vars)
    if (!info$uses_abun && length(info$cross_vars) > 0L) {
      .mgnet_abort_cross_without_abun(info$expr, info$cross_vars, verb, ctx)
    }
  }
  
  invisible(NULL)
}

# =============================================================================
# Long abundance builder (abun + optional cross vars)
# =============================================================================

.mgnet_long_abundance_with_cross <- function(object, ctx, abun_keys, cross_vars) {
  if (length(abun_keys) == 0L) return(tibble::tibble())
  
  # Long table for a single mgnet: sample_id x taxa_id + requested matrices
  .internal_long_abun <- function(obj, keys) {
    out <- tidyr::expand_grid(sample_id = sample_id(obj), taxa_id = taxa_id(obj))
    for (k in keys) {
      mat_long <- methods::slot(obj, k) |>
        tibble::as_tibble(rownames = "sample_id") |>
        tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = k)
      out <- dplyr::left_join(out, mat_long, by = c("sample_id", "taxa_id"))
    }
    out
  }
  
  # mgnet vs mgnets
  if (!ctx$is_list) {
    long_tbl <- .internal_long_abun(object, abun_keys)
  } else {
    long_tbl <- purrr::imap(
      object,
      ~ tibble::add_column(.internal_long_abun(.x, abun_keys), mgnet = .y, .before = 1)
    ) |>
      purrr::list_rbind()
  }
  
  # Attach cross vars via cross key (and mgnet if list)
  if (length(cross_vars) > 0L) {
    keep_cols <- unique(c(ctx$join_cross_key, cross_vars))
    cross_tbl <- ctx$get_cross_tbl(object) |>
      dplyr::select(tidyselect::all_of(keep_cols))
    long_tbl <- dplyr::left_join(long_tbl, cross_tbl, by = ctx$join_cross_key)
  }
  
  long_tbl
}

.mgnet_join_target_tbl <- function(long_tbl, target_tbl, ctx) {
  dplyr::left_join(long_tbl, target_tbl, by = ctx$join_target_key)
}

# =============================================================================
# Grouping semantics
# =============================================================================

.mgnet_groups <- function(object, ctx, is_abun_expr) {
  g <- get_group_mgnet(object)
  
  if (is.null(g)) {
    if (isTRUE(is_abun_expr)) {
      g <- if (!ctx$is_list) ctx$target_key else c("mgnet", ctx$target_key)
      return(rlang::syms(g))
    }
    return(rlang::syms(character(0)))
  }
  
  rlang::syms(g)
}

# =============================================================================
# Apply verbs (always returns ungrouped tbl)
# =============================================================================

.mgnet_apply_mutate <- function(tbl, expr_q, groups) {
  tbl |>
    dplyr::group_by(!!!groups) |>
    dplyr::mutate(!!!expr_q) |>
    dplyr::ungroup()
}

.mgnet_apply_filter <- function(tbl, exprs_q, groups) {
  tbl |>
    dplyr::group_by(!!!groups) |>
    dplyr::filter(!!!exprs_q) |>
    dplyr::ungroup()
}

# =============================================================================
# Cleanup after abun-aware evaluation (STRICT)
# - Drop helper cols
# - Require exactly 1 row per target id (else error)
# =============================================================================

.mgnet_cleanup_after_abun <- function(tbl, ctx,
                                      abun_drop  = c("abun", "rela", "norm"),
                                      cross_drop = character(0),
                                      verb = "mutate",
                                      expr_label = NULL) {
  
  out <- tbl |>
    dplyr::select(-tidyselect::any_of(c(ctx$cross_key, abun_drop, cross_drop))) |>
    dplyr::distinct()
  
  dup <- out |>
    dplyr::count(!!!rlang::syms(ctx$target_key), name = "n") |>
    dplyr::filter(.data$n > 1L)
  
  if (nrow(dup) > 0L) {
    show <- dup |>
      dplyr::select(!!!rlang::syms(ctx$target_key)) |>
      dplyr::slice_head(n = 8L)
    
    ex_keys <- paste(apply(show, 1, paste, collapse = ": "), collapse = ", ")
    
    cli::cli_abort(
      c(
        "x" = "Result is not unique per {ctx$target_key} after abundance-aware evaluation.",
        "i" = "You are calling {.fn {paste0(verb, '_', ctx$target)}} but the result still depends on the other dimension ({ctx$cross_key}).",
        "i" = "Collapse across that dimension (e.g. summaries like sum/mean/max).",
        "*" = if (!is.null(expr_label)) paste0("Expression: ", expr_label) else NULL,
        "*" = "Example offending keys: {.val {ex_keys}}"
      ),
      class = c(paste0("mgnet_", verb, "_", ctx$target, "_not_unique"), "mgnet_error")
    )
  }
  
  return(out)
}

.mgnet_assign_target_tbl <- function(object, ctx, tbl, verb, expr_label) {
  tryCatch(
    {
      if (ctx$target == "meta") meta(object) <- tbl else taxa(object) <- tbl
      object
    },
    error = function(err) .mgnet_abort_assign_failed(err, verb, ctx, expr_label)
  )
}

# =============================================================================
# Filter finalize (AND + reorder + subset) -- strict: NA -> error
# =============================================================================

.mgnet_filter_finalize <- function(object, ctx, filtered_lists, .ungroup = FALSE, expr_label = "AND across filters") {
  tryCatch(
    {
      key      <- ctx$target_key
      join_key <- ctx$join_target_key
      
      # Reorder by original ids; fail if any keep id is missing
      .reorder_strict <- function(keep, orig) {
        idx <- match(keep, orig)
        if (anyNA(idx)) {
          bad <- unique(keep[is.na(idx)])
          cli::cli_abort(
            c(
              "x" = "Filtering produced invalid {ctx$target} identifiers.",
              "*" = "Missing in original object: {.val {bad}}",
              "*" = "Stage: {.val {expr_label}}"
            ),
            class = c("mgnet_filter_invalid_ids", "mgnet_error")
          )
        }
        orig[idx]
      }
      
      # mgnet: vectors of ids -> intersect + reorder + subset
      if (!ctx$is_list) {
        keep <- purrr::reduce(filtered_lists, intersect)
        orig <- if (ctx$target == "meta") sample_id(object) else taxa_id(object)
        keep <- .reorder_strict(keep, orig)
        
        if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
        return(if (ctx$target == "meta") object[keep, ] else object[, keep])
      }
      
      # mgnets: tibbles with join_key -> semi_join AND + subset per mgnet
      final_tbl <- purrr::reduce(filtered_lists, dplyr::semi_join, by = join_key)
      
      for (nm in names(object)) {
        keep <- dplyr::filter(final_tbl, mgnet == nm) |>
          dplyr::pull(key)
        
        orig <- if (ctx$target == "meta") sample_id(object[[nm]]) else taxa_id(object[[nm]])
        keep <- .reorder_strict(keep, orig)
        
        object[[nm]] <- if (ctx$target == "meta") object[[nm]][keep, ] else object[[nm]][, keep]
      }
      
      validObject(object)
      if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
      object
    },
    error = function(err) .mgnet_abort_filter_failed(err, ctx, stage = "AND + reorder + subset", expr_label = expr_label)
  )
}

# =============================================================================
# Engine: MUTATE
# =============================================================================

.mgnet_mutate_engine <- function(object, ..., target = c("meta", "taxa"), .ungroup = FALSE) {
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  
  ctx <- .mgnet_ctx(object, target)
  .mgnet_abort_missing_slot(object, ctx)
  .mgnet_validate_exprs(object, ctx, exprs, verb = "mutate")
  
  # Build long table once (only if any expression uses abundance)
  info_all <- .mgnet_info_exprs(object, ctx, exprs)
  
  long_tbl <- if (isTRUE(info_all$uses_abun)) {
    .mgnet_long_abundance_with_cross(object, ctx,
                                     abun_keys  = info_all$abun_keys,
                                     cross_vars = info_all$cross_vars)
  } else {
    tibble::tibble()
  }
  
  for (i in seq_along(exprs)) {
    q          <- exprs[[i]]
    q_info     <- .mgnet_info_quo(object, ctx, q)
    expr_label <- rlang::as_label(q)
    
    target_tbl <- ctx$get_target_tbl(object)
    groups     <- .mgnet_groups(object, ctx, is_abun_expr = q_info$uses_abun)

    if (isTRUE(q_info$uses_abun)) {
      tbl <- .mgnet_join_target_tbl(long_tbl, target_tbl, ctx)
      tbl <- .mgnet_apply_mutate(tbl, exprs[i], groups)
      tbl <- .mgnet_cleanup_after_abun(
        tbl, ctx,
        cross_drop = info_all$cross_vars,
        verb = "mutate",
        expr_label = expr_label
      )
    } else {
      tbl <- .mgnet_apply_mutate(target_tbl, exprs[i], groups)
    }
    
    object <- .mgnet_assign_target_tbl(object, ctx, tbl, verb = "mutate", expr_label = expr_label)
  }
  
  if (isTRUE(.ungroup)) object <- ungroup_mgnet(object)
  object
}

# =============================================================================
# Engine: FILTER
# =============================================================================

.mgnet_filter_engine <- function(object, ..., target = c("meta", "taxa"), .ungroup = FALSE) {
  exprs <- rlang::enquos(...)
  if (length(exprs) == 0L) return(object)
  
  ctx <- .mgnet_ctx(object, target)
  .mgnet_abort_missing_slot(object, ctx)
  .mgnet_validate_exprs(object, ctx, exprs, verb = "filter")
  
  target_tbl <- ctx$get_target_tbl(object)
  
  # Build long table once (only if any expression uses abundance)
  info_all <- .mgnet_info_exprs(object, ctx, exprs)
  
  long_tbl <- if (isTRUE(info_all$uses_abun)) {
    .mgnet_long_abundance_with_cross(object, ctx,
                                     abun_keys  = info_all$abun_keys,
                                     cross_vars = info_all$cross_vars) |>
      .mgnet_join_target_tbl(target_tbl, ctx)
  } else {
    tibble::tibble()
  }
  
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    q          <- exprs[[i]]
    q_info     <- .mgnet_info_quo(object, ctx, q)
    expr_label <- rlang::as_label(q)
    
    groups <- .mgnet_groups(object, ctx, is_abun_expr = q_info$uses_abun)
    
    if (isTRUE(q_info$uses_abun)) {
      tbl <- .mgnet_apply_filter(long_tbl, exprs[i], groups)
      tbl <- .mgnet_cleanup_after_abun(
        tbl, ctx,
        cross_drop = info_all$cross_vars,
        verb = "filter",
        expr_label = expr_label
      )
    } else {
      tbl <- .mgnet_apply_filter(target_tbl, exprs[i], groups)
    }
    
    # Collect ids to AND later
    filtered_lists[[i]] <- if (!ctx$is_list) {
      dplyr::pull(tbl, ctx$target_key)
    } else {
      dplyr::select(tbl, tidyselect::all_of(ctx$join_target_key))
    }
  }
  
  .mgnet_filter_finalize(object, ctx, filtered_lists, .ungroup = .ungroup, expr_label = "AND across filters")
}

# =============================================================================
# Public API: generics + methods
# =============================================================================

#' @rdname mutate_mgnet
#' @export
setGeneric("mutate_meta", function(object, ..., .ungroup = FALSE) standardGeneric("mutate_meta"))

#' @rdname mutate_mgnet
#' @export
setMethod("mutate_meta", "mgnet", function(object, ..., .ungroup = FALSE) {
  .mgnet_mutate_engine(object, ..., target = "meta", .ungroup = .ungroup)
})

#' @rdname mutate_mgnet
#' @export
setMethod("mutate_meta", "mgnets", function(object, ..., .ungroup = FALSE) {
  .mgnet_mutate_engine(object, ..., target = "meta", .ungroup = .ungroup)
})

#' @rdname mutate_mgnet
#' @export
setGeneric("mutate_taxa", function(object, ..., .ungroup = FALSE) standardGeneric("mutate_taxa"))

#' @rdname mutate_mgnet
#' @export
setMethod("mutate_taxa", "mgnet", function(object, ..., .ungroup = FALSE) {
  .mgnet_mutate_engine(object, ..., target = "taxa", .ungroup = .ungroup)
})

#' @rdname mutate_mgnet
#' @export
setMethod("mutate_taxa", "mgnets", function(object, ..., .ungroup = FALSE) {
  .mgnet_mutate_engine(object, ..., target = "taxa", .ungroup = .ungroup)
})

#' @rdname filter_mgnet
#' @export
setGeneric("filter_meta", function(object, ..., .ungroup = FALSE) standardGeneric("filter_meta"))

#' @rdname filter_mgnet
#' @export
setMethod("filter_meta", "mgnet", function(object, ..., .ungroup = FALSE) {
  .mgnet_filter_engine(object, ..., target = "meta", .ungroup = .ungroup)
})

#' @rdname filter_mgnet
#' @export
setMethod("filter_meta", "mgnets", function(object, ..., .ungroup = FALSE) {
  .mgnet_filter_engine(object, ..., target = "meta", .ungroup = .ungroup)
})

#' @rdname filter_mgnet
#' @export
setGeneric("filter_taxa", function(object, ..., .ungroup = FALSE) standardGeneric("filter_taxa"))

#' @rdname filter_mgnet
#' @export
setMethod("filter_taxa", "mgnet", function(object, ..., .ungroup = FALSE) {
  .mgnet_filter_engine(object, ..., target = "taxa", .ungroup = .ungroup)
})

#' @rdname filter_mgnet
#' @export
setMethod("filter_taxa", "mgnets", function(object, ..., .ungroup = FALSE) {
  .mgnet_filter_engine(object, ..., target = "taxa", .ungroup = .ungroup)
})
