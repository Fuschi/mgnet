# Internal function to check the assigned elements for the setter methods of mgnets
#------------------------------------------------------------------------------#

#' Check list in assign methods for mgnets
#'
#' @description
#' Internal validator ensuring that a named list supplied to an `mgnets`
#' setter has the correct length and name set (order may differ).
#'
#' @param object An `mgnets` object.
#' @param value  A named list intended to be assigned to the `mgnets`.
#'
#' @return `invisible(TRUE)` if all checks pass; otherwise throws an error.
#' @keywords internal
#' @export
is_list_mgnets_assign <- function(object, value) {
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  # 1) Type check on object ----------------------------------------------------
  if (!inherits(object, "mgnets")) {
    cli::cli_abort("{.arg {objectName}} must be an {.cls mgnets} object.")
  }
  
  # 2) Basic type/length checks on value ---------------------------------------
  if (!is.list(value)) {
    cli::cli_abort("{.arg {valueName}} must be a list.")
  }
  
  obj_len <- length(object@mgnets)
  val_len <- length(value)
  if (obj_len != val_len) {
    cli::cli_abort(
      "Lengths must match: {obj_len} in {.arg {objectName}}@mgnets vs {val_len} in {.arg {valueName}}."
    )
  }
  
  # 3) Names: must exist, non-empty, no duplicates -----------------------------
  val_names <- names(value)
  if (is.null(val_names) || any(is.na(val_names)) || any(val_names == "")) {
    cli::cli_abort("{.arg {valueName}} must be a named list with non-empty names.")
  }
  
  dup_idx <- anyDuplicated(val_names)
  if (dup_idx) {
    dups <- unique(val_names[duplicated(val_names)])
    dups_fmt <- paste0("{", paste(dups, collapse = "}, {"), "}")
    cli::cli_abort("{.arg {valueName}} cannot contain duplicated names: {dups_fmt}.")
  }
  
  # 4) Name set equality (order can differ) ------------------------------------
  obj_names <- names(object@mgnets)
  if (!is.null(obj_names)) {
    miss_in_value  <- setdiff(obj_names, val_names)
    extra_in_value <- setdiff(val_names, obj_names)
    
    if (length(miss_in_value) || length(extra_in_value)) {
      miss_fmt  <- if (length(miss_in_value)) paste(miss_in_value, collapse = ", ") else "none"
      extra_fmt <- if (length(extra_in_value)) paste(extra_in_value, collapse = ", ") else "none"
      cli::cli_abort(c(
        "Names of {.arg {valueName}} must match names of {.arg {objectName}}@mgnets (order can differ).",
        "x" = "Missing in {.arg {valueName}}: {miss_fmt}",
        "x" = "Extra in {.arg {valueName}}: {extra_fmt}"
      ))
    }
  }
  
  invisible(TRUE)
}

#' Check Tibble in Assign Methods for mgnets
#'
#' @description
#' Internal validator for a tibble/data frame replacing all sample/taxa
#' metadata across an `mgnets` object. Requires `mgnet` and `sample_id`/`taxa_id`
#' columns and forbids row names.
#'
#' @param object An `mgnets` object.
#' @param value  A tibble/data frame with columns `mgnet` and either
#'   `sample_id` (if `sample_or_taxa = "sample"`) or `taxa_id` (if `"taxa"`).
#' @param sample_or_taxa Character, `"sample"` or `"taxa"`.
#'
#' @return `invisible(TRUE)` if all checks pass; otherwise throws an error.
#' @keywords internal
#' @importFrom cli cli_abort
#' @importFrom tibble has_rownames
#' @export
is_assign_mgnets_tbl <- function(object, value, sample_or_taxa) {
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  sample_or_taxa <- match.arg(sample_or_taxa, c("sample", "taxa"))
  
  # 1) Object type -------------------------------------------------------------
  if (!inherits(object, "mgnets")) {
    cli::cli_abort("{.arg {objectName}} must be an {.cls mgnets} object.")
  }
  
  # 2) Value type --------------------------------------------------------------
  if (!(inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame.")
  }
  
  # 3) Required columns --------------------------------------------------------
  required_cols <- if (sample_or_taxa == "sample") c("mgnet", "sample_id") else c("mgnet", "taxa_id")
  value_cols <- colnames(value)
  missing_cols <- setdiff(required_cols, value_cols)
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Missing required columns in {.arg {valueName}}.",
      "x" = "Required: {paste(required_cols, collapse = ', ')}",
      "x" = "Missing: {paste(missing_cols, collapse = ', ')}"
    ))
  }
  
  # 4) Rownames must NOT be set ------------------------------------------------
  if (tibble::has_rownames(value)) {
    cli::cli_abort("Row names must not be set in {.arg {valueName}}.")
  }
  
  # 5) 'mgnet' values must match names(object) ----------------------------------
  unique_mgnets  <- unique(value$mgnet)
  object_mgnets  <- names(object)
  missing_mgnets <- setdiff(unique_mgnets, object_mgnets)
  if (length(missing_mgnets) > 0L) {
    cli::cli_abort(c(
      "Some {.val mgnet} values in {.arg {valueName}} are not present in {.arg {objectName}}.",
      "x" = "Not found: {paste(missing_mgnets, collapse = ', ')}",
      "i" = "Allowed values: {paste(object_mgnets, collapse = ', ')}"
    ))
  }
  
  invisible(TRUE)
}


#' Check Tibble in Assign Methods for mgnets links
#'
#' @description
#' Internal validator for a tibble/data frame replacing all links 
#' metadata across an `mgnets` object. Requires `mgnet`,`from`,`to` and `link_id`
#' columns and forbids row names.
#'
#' @param object An `mgnets` object.
#' @param value  A tibble/data frame with columns `mgnet`,`from`,`to` and `link_id` 
#' and the additional information
#'
#' @return `invisible(TRUE)` if all checks pass; otherwise throws an error.
#' @keywords internal
#' @export
is_assign_mgnets_link_tbl <- function(object, value) {
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))

  # 1) Object type -------------------------------------------------------------
  if (!inherits(object, "mgnets")) {
    cli::cli_abort("{.arg {objectName}} must be an {.cls mgnets} object.")
  }
  
  # 2) Value type --------------------------------------------------------------
  if (!(inherits(value, "data.frame") || inherits(value, "tbl_df"))) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame.")
  }
  
  # 3) Required columns --------------------------------------------------------
  required_cols <- c("mgnet","from","to","link_id")
  value_cols <- colnames(value)
  missing_cols <- setdiff(required_cols, value_cols)
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Missing required columns in {.arg {valueName}}.",
      "x" = "Required: {paste(required_cols, collapse = ', ')}",
      "x" = "Missing: {paste(missing_cols, collapse = ', ')}"
    ))
  }
  
  # 4) Rownames must NOT be set ------------------------------------------------
  if (tibble::has_rownames(value)) {
    cli::cli_abort("Row names must not be set in {.arg {valueName}}.")
  }
  
  # 5) 'mgnet' values must match names(object) ----------------------------------
  unique_mgnets  <- unique(value$mgnet)
  object_mgnets  <- names(object)
  missing_mgnets <- setdiff(unique_mgnets, object_mgnets)
  if (length(missing_mgnets) > 0L) {
    cli::cli_abort(c(
      "Some {.val mgnet} values in {.arg {valueName}} are not present in {.arg {objectName}}.",
      "x" = "Not found: {paste(missing_mgnets, collapse = ', ')}",
      "i" = "Allowed values: {paste(object_mgnets, collapse = ', ')}"
    ))
  }
  
  invisible(TRUE)
}