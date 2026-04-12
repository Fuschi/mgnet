#' @include class-mgnet.R class-mgnets.R
NULL

# Internal utility ---------------------------------------------------------

.build_split_name <- function(group_df) {
  paste(
    vapply(
      group_df[1, , drop = FALSE],
      function(x) as.character(x[[1]]),
      character(1)
    ),
    collapse = "-"
  )
}

.split_mgnet_by_tbl <- function(object, split_tbl, full_tbl, margin = c("sample", "taxa")) {
  
  margin <- match.arg(margin)
  
  groups <- split_tbl %>%
    dplyr::distinct()
  
  mgnet_list <- vector("list", nrow(groups))
  group_names <- character(nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    group_df <- groups[i, , drop = FALSE]
    
    filtered_tbl <- full_tbl %>%
      dplyr::semi_join(group_df, by = colnames(group_df))
    
    subsetted_mgnet <- switch(
      margin,
      sample = object[filtered_tbl$sample_id, , drop = FALSE],
      taxa   = object[, filtered_tbl$taxa_id, drop = FALSE]
    )
    
    group_names[i] <- .build_split_name(group_df)
    mgnet_list[[i]] <- subsetted_mgnet
  }
  
  if (anyDuplicated(group_names) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Generated group names are not unique.",
        "i" = "Try selecting splitting variables whose value combinations produce unique names."
      )
    )
  }
  
  names(mgnet_list) <- group_names
  mgnets(mgnet_list)
}

#' Split an `mgnet` object into multiple `mgnet` objects
#'
#' @description
#' Split an `mgnet` object into multiple `mgnet` objects according to values of
#' selected metadata columns. The resulting objects are returned inside an
#' `mgnets` container.
#'
#' `split_meta()` performs the split using **sample-level metadata** returned
#' by [meta()], while `split_taxa()` performs the split using **taxa-level
#' metadata** returned by [taxa()].
#'
#' @param object An object of class `mgnet`.
#' @param ... One or more metadata columns used to define the split.
#'
#' @details
#' Each unique combination of the selected columns defines a group. A new
#' `mgnet` object is created for each group.
#'
#' The names of the resulting objects are constructed by concatenating the
#' values of the selected columns using `"-"`.
#'
#' Group names must be unique.
#'
#' @return
#' An object of class `mgnets` containing one `mgnet` per group.
#'
#' @examples
#' \dontrun{
#' split_meta(x, Country)
#' split_meta(x, Country, Surface)
#'
#' split_taxa(x, Phylum)
#' split_taxa(x, Phylum, Genus)
#' }
#'
#' @name split_mgnet
NULL

#' @rdname split_mgnet
#' @export
setGeneric("split_meta", function(object, ...) standardGeneric("split_meta"))

#' @rdname split_mgnet
#' @export
setMethod("split_meta", "mgnet", function(object, ...) {
  
  if (miss_sample(object)) {
    cli::cli_abort("No sample available.")
  }
  
  split_cols <- rlang::enquos(...)
  
  if (length(split_cols) == 0L) {
    cli::cli_abort("No columns specified for splitting. Please specify one or more columns.")
  }
  
  meta_tbl <- meta(object, .fmt = "tbl")
  
  split_tbl <- meta_tbl %>%
    dplyr::select(!!!split_cols)
  
  .split_mgnet_by_tbl(
    object   = object,
    split_tbl = split_tbl,
    full_tbl  = meta_tbl,
    margin    = "sample"
  )
})

#' @rdname split_mgnet
#' @export
setGeneric("split_taxa", function(object, ...) standardGeneric("split_taxa"))

#' @rdname split_mgnet
#' @export
setMethod("split_taxa", "mgnet", function(object, ...) {
  
  if (miss_taxa(object)) {
    cli::cli_abort("No taxa available.")
  }
  
  split_cols <- rlang::enquos(...)
  
  if (length(split_cols) == 0L) {
    cli::cli_abort("No columns specified for splitting. Please specify one or more columns.")
  }
  
  taxa_tbl <- taxa(object, .fmt = "tbl")
  
  split_tbl <- taxa_tbl %>%
    dplyr::select(!!!split_cols)
  
  .split_mgnet_by_tbl(
    object    = object,
    split_tbl = split_tbl,
    full_tbl  = taxa_tbl,
    margin    = "taxa"
  )
})