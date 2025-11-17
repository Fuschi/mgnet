#' @include class-mgnet.R class-mgnets.R class-getters.R
NULL

# -------------------------------------------------------------------------
# INTERNAL HELPERS (not exported)
# -------------------------------------------------------------------------

#' @noRd
.match_has_miss_fmt <- function(x) {
  allowed <- c("list", "any", "all")
  if (is.null(x)) return("list")
  x <- tolower(as.character(x))
  if (!x %in% allowed) {
    cli::cli_abort(c(
      "x" = "Invalid {.arg .fmt}: {.val {x}}.",
      "v" = "Allowed values are {.val 'list'}, {.val 'any'}, {.val 'all'}."
    ))
  }
  x
}

#' @noRd
.is_nonempty <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(nrow(x) > 0 || ncol(x) > 0)
  }
  length(x) > 0
}

#' @noRd
.has_slot_mgnet <- function(object, slot) {
  if (!methods::is(object, "mgnet")) {
    cli::cli_abort("'.has_slot_mgnet' expects an 'mgnet' object.")
  }
  .is_nonempty(methods::slot(object, slot))
}

#' @noRd
.has_slot_mgnets <- function(object, slot) {
  if (!methods::is(object, "mgnets")) {
    cli::cli_abort("'.has_slot_mgnets' expects an 'mgnets' object.")
  }
  sapply(object@mgnets, function(x) .has_slot_mgnet(x, slot),
         simplify = TRUE, USE.NAMES = TRUE)
}

#' @noRd
.slot_exists <- function(object, slot_name) {
  slot_name %in% methods::slotNames(object)
}

#' @importFrom methods slot slotNames
NULL


# =============================================================================
# NSAMPLE
# =============================================================================

#' Get Number of Samples
#'
#' Returns the number of samples in an `mgnet` or `mgnets` object.
#'
#' For an empty object, returns `0L`. For `mgnets`, returns a named integer
#' vector with one value per element.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @return Integer (for `mgnet`) or named integer vector (for `mgnets`).
#' @export
#' @name nsample
#' @aliases nsample,mgnet-method nsample,mgnets-method
setGeneric("nsample", function(object) standardGeneric("nsample"))

#' @export
setMethod("nsample", "mgnet", function(object) {
  if (.has_slot_mgnet(object, "abun")) return(nrow(object@abun))
  if (.has_slot_mgnet(object, "rela")) return(nrow(object@rela))
  if (.has_slot_mgnet(object, "norm")) return(nrow(object@norm))
  if (.has_slot_mgnet(object, "meta")) return(nrow(object@meta))
  0L
})

#' @export
setMethod("nsample", "mgnets", function(object) {
  sapply(object@mgnets, nsample, simplify = TRUE, USE.NAMES = TRUE)
})


# =============================================================================
# NTAXA
# =============================================================================

#' Get Number of Taxa
#'
#' Returns the number of taxa in an `mgnet` or `mgnets` object.
#'
#' For an empty object, returns `0L`. For `mgnets`, returns a named integer
#' vector with one value per element.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @return Integer (for `mgnet`) or named integer vector (for `mgnets`).
#' @export
#' @name ntaxa
#' @aliases ntaxa,mgnet-method ntaxa,mgnets-method
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))

#' @export
setMethod("ntaxa", "mgnet", function(object) {
  if (.has_slot_mgnet(object, "abun")) return(ncol(object@abun))
  if (.has_slot_mgnet(object, "taxa")) return(nrow(object@taxa))
  if (.has_slot_mgnet(object, "rela")) return(ncol(object@rela))
  if (.has_slot_mgnet(object, "norm")) return(ncol(object@norm))
  0L
})

#' @export
setMethod("ntaxa", "mgnets", function(object) {
  sapply(object@mgnets, ntaxa, simplify = TRUE, USE.NAMES = TRUE)
})


# =============================================================================
# SAMPLE_ID
# =============================================================================

#' Get Sample IDs
#'
#' Extracts sample IDs from an `mgnet` or `mgnets` object. Returns an empty
#' character vector when no samples are available.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @return Character vector (for `mgnet`), or a named list of character vectors (for `mgnets`).
#' @export
#' @name sample_id
#' @aliases sample_id,mgnet-method sample_id,mgnets-method
setGeneric("sample_id", function(object) standardGeneric("sample_id"))

#' @export
setMethod("sample_id", "mgnet", function(object) {
  if (.has_slot_mgnet(object, "abun")) return(rownames(object@abun))
  if (.has_slot_mgnet(object, "meta")) return(rownames(object@meta))
  if (.has_slot_mgnet(object, "norm")) return(rownames(object@norm))
  if (.has_slot_mgnet(object, "rela")) return(rownames(object@rela))
  character(0)
})

#' @export
setMethod("sample_id", "mgnets", function(object) {
  sapply(object@mgnets, sample_id, simplify = FALSE, USE.NAMES = TRUE)
})


# =============================================================================
# TAXA_ID
# =============================================================================

#' Get Taxa IDs
#'
#' Extracts taxa IDs from an `mgnet` or `mgnets` object. Returns an empty
#' character vector when no taxa are available.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @return Character vector (for `mgnet`), or a named list of character vectors (for `mgnets`).
#' @export
#' @name taxa_id
#' @aliases taxa_id,mgnet-method taxa_id,mgnets-method
setGeneric("taxa_id", function(object) standardGeneric("taxa_id"))

#' @export
setMethod("taxa_id", "mgnet", function(object) {
  if (.has_slot_mgnet(object, "abun")) return(colnames(object@abun))
  if (.has_slot_mgnet(object, "taxa")) return(rownames(object@taxa))
  if (.has_slot_mgnet(object, "norm")) return(colnames(object@norm))
  if (.has_slot_mgnet(object, "rela")) return(colnames(object@rela))
  character(0)
})

#' @export
setMethod("taxa_id", "mgnets", function(object) {
  sapply(object@mgnets, taxa_id, simplify = FALSE, USE.NAMES = TRUE)
})


# =============================================================================
# NCOMM (number of communities) + COMM_ID
# =============================================================================

#' Get Number of Communities
#'
#' Returns the number of detected communities in an `mgnet` or each element of an `mgnets`.
#' Returns `0L` when communities are absent.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @return Integer (for `mgnet`) or named integer vector (for `mgnets`).
#' @export
#' @name ncomm
#' @aliases ncomm,mgnet-method ncomm,mgnets-method
#' @importFrom igraph sizes
setGeneric("ncomm", function(object) standardGeneric("ncomm"))

#' @export
setMethod("ncomm", "mgnet", function(object) {
  if (.has_slot_mgnet(object, "comm")) {
    return(length(igraph::sizes(object@comm)))
  }
  0L
})

#' @export
setMethod("ncomm", "mgnets", function(object) {
  sapply(object@mgnets, ncomm, simplify = TRUE, USE.NAMES = TRUE)
})

#' Retrieve Community IDs
#'
#' Returns community membership for taxa in an `mgnet` object or for each element of an `mgnets`.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @param .fmt Output format. One of:
#' \itemize{
#'   \item `"list"`: a named character vector (for `mgnet`) or a named list (for `mgnets`) of community IDs.
#'   \item `"tbl"`: a tibble with columns `taxa_id`, `comm_id`, and `mgnet` (the latter only for `mgnets`).
#' }
#' Default is `"list"`.
#'
#' @return Character vector, list, or tibble depending on `.fmt` and the class of `object`.
#' Returns empty structures when communities are absent.
#'
#' @importFrom igraph membership
#' @importFrom tibble tibble
#' @importFrom purrr imap list_rbind
#' @importFrom stats setNames
#' @export
#' @name comm_id
#' @aliases comm_id,mgnet-method comm_id,mgnets-method
setGeneric("comm_id", function(object, .fmt = "list") standardGeneric("comm_id"))

#' @export
setMethod("comm_id", "mgnet", function(object, .fmt = "list") {
  
  if (!.fmt %in% c("list", "tbl")) {
    cli::cli_abort(c(
      "x" = "Invalid {.arg .fmt} argument: {.val {.fmt}}.",
      "v" = "Allowed values are {.val 'list'} or {.val 'tbl'}."
    ))
  }
  
  if (.has_slot_mgnet(object, "comm")) {
    memb <- as.character(igraph::membership(object@comm))
    ids  <- taxa_id(object)
    if (.fmt == "list") {
      stats::setNames(memb, ids)
    } else {
      tibble::tibble(taxa_id = ids, comm_id = memb)
    }
  } else {
    if (.fmt == "tbl") {
      tibble::tibble(taxa_id = character(0), comm_id = character(0))
    } else {
      stats::setNames(character(0), character(0))
    }
  }
})

#' @export
setMethod("comm_id", "mgnets", function(object, .fmt = "list") {
  
  if (!.fmt %in% c("list", "tbl")) {
    cli::cli_abort(c(
      "x" = "Invalid {.arg .fmt} argument: {.val {.fmt}}.",
      "v" = "Allowed values are {.val 'list'} or {.val 'tbl'}."
    ))
  }
  
  if (.fmt == "list") {
    lapply(object@mgnets, function(x) comm_id(x, .fmt = "list"))
  } else {
    object@mgnets |>
      purrr::imap(function(mgnet_obj, name) {
        if (.has_slot_mgnet(mgnet_obj, "comm")) {
          tibble::tibble(
            mgnet    = name,
            taxa_id = taxa_id(mgnet_obj),
            comm_id = as.character(igraph::membership(mgnet_obj@comm))
          )
        } else {
          tibble::tibble(
            mgnet    = name,
            taxa_id = character(0),
            comm_id = character(0)
          )
        }
      }) |>
      purrr::list_rbind()
  }
})


# =============================================================================
# LINK ID
# =============================================================================

#' Retrieve Link IDs
#'
#' @description
#' Return the `link_id` edge attribute from the network stored in an `mgnet`
#' object, or a named list of `link_id` vectors for each element of an `mgnets`.
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @param selected Logical. If `TRUE` (default), return link IDs from
#'   the currently **selected** subgraph (if any); if `FALSE`, return link IDs
#'   from the full network.
#'
#' @return
#' - For an `mgnet`: an atmgnet vector (usually `character`) of link IDs taken
#'   from the edge attribute `link_id` of `netw(object, selected = selected)`.
#' - For an `mgnets`: a **named list** with one element per contained `mgnet`,
#'   each element being the vector described above.
#'
#' @details
#' This function is a light wrapper around `igraph::edge_attr()` applied to the
#' network returned by `netw(object, selected = selected)`. It assumes that the
#' network carries an edge attribute named `link_id`.
#'
#' @export
#' @name link_id
#' @aliases link_id,mgnet-method link_id,mgnets-method
setGeneric("link_id", function(object, selected = TRUE) standardGeneric("link_id"))

#' @export
setMethod("link_id", "mgnet", function(object, selected = TRUE) {

  if (!(is.logical(selected) && length(selected) == 1L && !is.na(selected))) {
    cli::cli_abort(c(
      "x" = "{.arg selected} must be a single {.cls logical} (TRUE/FALSE) and not NA.",
      "i" = "Got class {.cls {class(selected)[1]}} with length {length(selected)}."
    ))
  }
  
  if(isTRUE(selected) && are_selected_links(object)){
    
    link_id <- igraph::edge_attr(object@netw, "link_id")
    return(link_id[link_id %in% get_selected_links(object)])
    
  } else {
    
    return(igraph::edge_attr(object@netw, "link_id"))
    
  }
  
})

#' @export
setMethod("link_id", "mgnets", function(object, selected = TRUE) {
  
  if (!(is.logical(selected) && length(selected) == 1L && !is.na(selected))) {
    cli::cli_abort(c(
      "x" = "{.arg selected} must be a single {.cls logical} (TRUE/FALSE) and not NA.",
      "i" = "Got class {.cls {class(selected)[1]}} with length {length(selected)}."
    ))
  }
  
  sapply(object@mgnets, \(x) link_id(x, selected = selected),
         simplify = FALSE, USE.NAMES = TRUE)
  
})



# =============================================================================
# HAS / MISS METHODS
# =============================================================================

#' Presence/absence helpers for `mgnet` / `mgnets`
#'
#' A unified set of helpers to check whether samples, taxa, or specific slots
#' are present (`has_*`) or missing (`miss_*`) in `mgnet` and `mgnets` objects.
#'
#' @section Output format for `mgnets`:
#' The `.fmt` argument controls the output when `object` is an `mgnets`:
#' \itemize{
#'   \item `"list"` (default): named logical vector (one value per element)
#'   \item `"any"`: single logical, `TRUE` if at least one element matches
#'   \item `"all"`: single logical, `TRUE` if all elements match
#' }
#' The value is validated internally by `.match_has_miss_fmt()`.
#'
#' @section Semantics:
#' \tabular{ll}{
#'   `has_sample` / `miss_sample` \tab use `nsample(object) > 0` / `== 0` \cr
#'   `has_taxa`   / `miss_taxa`   \tab use `ntaxa(object) > 0` / `== 0` \cr
#'   `has_slot`   / `miss_slot`   \tab check that a slot exists and `length(slot) > 0` / `== 0` \cr
#'   `has_abun`   / `miss_abun`   \tab wrappers for `"abun"` \cr
#'   `has_rela`   / `miss_rela`   \tab wrappers for `"rela"` \cr
#'   `has_norm`   / `miss_norm`   \tab wrappers for `"norm"` \cr
#'   `has_meta`   / `miss_meta`   \tab wrappers for `"meta"` \cr
#'   `has_netw`   / `miss_netw`   \tab wrappers for `"network"` \cr
#'   `has_comm`   / `miss_comm`   \tab wrappers for `"comm"` \cr
#'   `has_metataxa` / `miss_metataxa` \tab presence/absence of at least one of `"taxa"` or `"comm"` \cr
#' }
#'
#' @section Edge cases:
#' \itemize{
#'   \item For an empty `mgnets` (length 0), `has_*` return `FALSE` and `miss_*` return `TRUE`.
#'   \item `has_slot`/`miss_slot` throw an error if the slot does not exist in class `mgnet`.
#' }
#'
#' @param object An object of class `mgnet` or `mgnets`.
#' @param slot_name (only for `has_slot`/`miss_slot`) Name of the slot to check.
#' @param .fmt (only for `mgnets`) One of `"list"`, `"any"`, `"all"`.
#' @return For `mgnet`: single logical. For `mgnets`: named logical vector (`"list"`)
#' or single logical (`"any"`, `"all"`).
#' @family mgnet-introspection
#' @seealso [nsample()], [ntaxa()]
#' @name has_miss-methods
#' @aliases
#' has_sample miss_sample
#' has_taxa miss_taxa
#' has_slot miss_slot
#' has_abun miss_abun
#' has_rela miss_rela
#' has_norm miss_norm
#' has_meta miss_meta
#' has_netw miss_netw
#' has_comm miss_comm
#' has_metataxa miss_metataxa
NULL

# --- sample -------------------------------------------------------------------

#' @rdname has_miss-methods
setGeneric("has_sample", function(object, .fmt = "list") standardGeneric("has_sample"))

setMethod("has_sample", "mgnet", function(object) {
  nsample(object) > 0
})

setMethod("has_sample", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_sample, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_sample", function(object, .fmt = "list") standardGeneric("miss_sample"))

setMethod("miss_sample", "mgnet", function(object) {
  nsample(object) == 0
})

setMethod("miss_sample", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_sample, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

# --- taxa ---------------------------------------------------------------------

#' @rdname has_miss-methods
setGeneric("has_taxa", function(object, .fmt = "list") standardGeneric("has_taxa"))

setMethod("has_taxa", "mgnet", function(object) {
  ntaxa(object) > 0
})

setMethod("has_taxa", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_taxa, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_taxa", function(object, .fmt = "list") standardGeneric("miss_taxa"))

setMethod("miss_taxa", "mgnet", function(object) {
  ntaxa(object) == 0
})

setMethod("miss_taxa", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_taxa, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

# --- generic slot -------------------------------------------------------------

#' @rdname has_miss-methods
setGeneric("has_slot", function(object, slot_name, .fmt = "list") standardGeneric("has_slot"))

setMethod("has_slot", "mgnet", function(object, slot_name, .fmt = "list") {
  if (!.slot_exists(object, slot_name)) {
    cli::cli_abort("Slot {.val {slot_name}} does not exist in {.cls mgnet}.")
  }
  length(methods::slot(object, slot_name)) > 0
})

setMethod("has_slot", "mgnets", function(object, slot_name, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_slot, slot_name = slot_name, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})


#' @rdname has_miss-methods
setGeneric("miss_slot", function(object, slot_name, .fmt = "list") standardGeneric("miss_slot"))

setMethod("miss_slot", "mgnet", function(object, slot_name, .fmt = "list") {
  if (!.slot_exists(object, slot_name)) {
    cli::cli_abort("Slot {.val {slot_name}} does not exist in {.cls mgnet}.")
  }
  length(methods::slot(object, slot_name)) == 0
})

setMethod("miss_slot", "mgnets", function(object, slot_name, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_slot, slot_name = slot_name, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

# --- wrappers for common slots ------------------------------------------------

#' @rdname has_miss-methods
setGeneric("has_abun", function(object, .fmt = "list") standardGeneric("has_abun"))
setMethod("has_abun", "mgnet",  function(object, .fmt = "list") has_slot(object, "abun"))
setMethod("has_abun", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_abun, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_abun", function(object, .fmt = "list") standardGeneric("miss_abun"))
setMethod("miss_abun", "mgnet",  function(object, .fmt = "list") miss_slot(object, "abun"))
setMethod("miss_abun", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_abun, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("has_rela", function(object, .fmt = "list") standardGeneric("has_rela"))
setMethod("has_rela", "mgnet",  function(object, .fmt = "list") has_slot(object, "rela"))
setMethod("has_rela", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_rela, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_rela", function(object, .fmt = "list") standardGeneric("miss_rela"))
setMethod("miss_rela", "mgnet",  function(object, .fmt = "list") miss_slot(object, "rela"))
setMethod("miss_rela", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_rela, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("has_norm", function(object, .fmt = "list") standardGeneric("has_norm"))
setMethod("has_norm", "mgnet",  function(object, .fmt = "list") has_slot(object, "norm"))
setMethod("has_norm", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_norm, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_norm", function(object, .fmt = "list") standardGeneric("miss_norm"))
setMethod("miss_norm", "mgnet",  function(object, .fmt = "list") miss_slot(object, "norm"))
setMethod("miss_norm", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_norm, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("has_meta", function(object, .fmt = "list") standardGeneric("has_meta"))
setMethod("has_meta", "mgnet",  function(object, .fmt = "list") has_slot(object, "meta"))
setMethod("has_meta", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_meta, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_meta", function(object, .fmt = "list") standardGeneric("miss_meta"))
setMethod("miss_meta", "mgnet",  function(object, .fmt = "list") miss_slot(object, "meta"))
setMethod("miss_meta", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_meta, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("has_netw", function(object, .fmt = "list") standardGeneric("has_netw"))
setMethod("has_netw", "mgnet",  function(object, .fmt = "list") has_slot(object, "network"))
setMethod("has_netw", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_netw, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_netw", function(object, .fmt = "list") standardGeneric("miss_netw"))
setMethod("miss_netw", "mgnet",  function(object, .fmt = "list") miss_slot(object, "netw"))
setMethod("miss_netw", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_netw, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("has_comm", function(object, .fmt = "list") standardGeneric("has_comm"))
setMethod("has_comm", "mgnet",  function(object, .fmt = "list") has_slot(object, "comm"))
setMethod("has_comm", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_comm, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_comm", function(object, .fmt = "list") standardGeneric("miss_comm"))
setMethod("miss_comm", "mgnet",  function(object, .fmt = "list") miss_slot(object, "comm"))
setMethod("miss_comm", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_comm, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

# --- composed helper: taxa OR comm -------------------------------------------

#' @rdname has_miss-methods
setGeneric("has_metataxa", function(object, .fmt = "list") standardGeneric("has_metataxa"))

setMethod("has_metataxa", "mgnet", function(object) {
  has_slot(object, "taxa") || has_slot(object, "comm")
})

setMethod("has_metataxa", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(FALSE)
  result <- sapply(object@mgnets, has_metataxa, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#' @rdname has_miss-methods
setGeneric("miss_metataxa", function(object, .fmt = "list") standardGeneric("miss_metataxa"))

setMethod("miss_metataxa", "mgnet", function(object) {
  miss_slot(object, "taxa") && miss_slot(object, "comm")
})

setMethod("miss_metataxa", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_has_miss_fmt(.fmt)
  if (length(object@mgnets) == 0) return(TRUE)
  result <- sapply(object@mgnets, miss_metataxa, USE.NAMES = TRUE)
  switch(.fmt, list = result, any = any(result), all = all(result))
})

#===============================================================================#
# META / TAXA / LINK / mgnet VARS
#===============================================================================#

#' @importFrom igraph edge_attr_names
NULL

#' @noRd
.match_vars_fmt <- function(x) {
  allowed <- c("list", "unique")
  if (is.null(x)) return("list")
  x <- tolower(as.character(x))
  if (!x %in% allowed) {
    cli::cli_abort(c(
      "x" = "Invalid {.arg .fmt}: {.val {x}}.",
      "v" = "Allowed values are {.val 'list'} or {.val 'unique'}."
    ))
  }
  x
}

#' Variable-name helpers for `mgnet` / `mgnets`
#'
#' Helpers to retrieve semantic variable names from an `mgnet` or an `mgnets`
#' collection:
#' \itemize{
#'   \item \code{meta_vars()}: sample-level variables. Returns \code{"sample_id"} when samples exist but \code{meta} is empty.
#'   \item \code{taxa_vars()}: taxa-level variables. Always includes \code{"taxa_id"} when taxa exist; adds \code{"comm_id"} if communities are present.
#'   \item \code{meta_taxa_vars()}: union of \code{meta_vars()} and \code{taxa_vars()}.
#'   \item \code{link_vars()}: edge-level variables from the network; returns \code{"link_id"} if a network exists but no edge attributes.
#'   \item \code{mgnet_vars()}: union of \code{meta_vars()}, \code{taxa_vars()}, and \code{link_vars()}.
#' }
#'
#' @section Output format for `mgnets`:
#' The \code{.fmt} argument controls the output when \code{object} is an \code{mgnets}:
#' \itemize{
#'   \item \code{"list"} (default): a named list, one character vector per element;
#'   \item \code{"unique"}: a single unique character vector across all elements.
#' }
#' In \code{"unique"} mode, \code{mgnet_vars()} and \code{meta_taxa_vars()} also include \code{"mgnet"} for convenience when stacking data.
#'
#' @section Missing data semantics:
#' \itemize{
#'   \item If no samples/taxa/network are present, the corresponding helper returns \code{character(0)}.
#'   \item If samples exist but \code{meta} is missing/empty, \code{meta_vars()} returns \code{"sample_id"}.
#'   \item If taxa exist but \code{taxa} is missing/empty, \code{taxa_vars()} returns \code{"taxa_id"} (and \code{"comm_id"} if communities exist).
#' }
#'
#' @param object An \code{mgnet} or \code{mgnets} object.
#' @param .fmt For \code{mgnets} only: one of \code{"list"} or \code{"unique"}.
#'
#' @return
#' \itemize{
#'   \item For \code{mgnet}: a character vector (possibly empty).
#'   \item For \code{mgnets}: named list of character vectors (\code{.fmt = "list"}) or a unique character vector (\code{.fmt = "unique"}).
#' }
#'
#' @family vars-helpers
#' @seealso \code{\link{has_meta}}, \code{\link{has_taxa}}, \code{\link{has_netw}}, \code{\link{has_comm}}
#' @name vars-helpers
#' @aliases
#' meta_vars taxa_vars meta_taxa_vars link_vars mgnet_vars
NULL

#------------------------------------------------------------------------------#
# meta_vars
#------------------------------------------------------------------------------#

#' @rdname vars-helpers
#' @export
setGeneric("meta_vars", function(object, .fmt = "list") standardGeneric("meta_vars"))

#' @rdname vars-helpers
#' @export
setMethod("meta_vars", "mgnet", function(object, .fmt = "list") {
  if (!has_sample(object)) return(character(0))
  if (!has_meta(object))   return("sample_id")
  unique(c("sample_id", colnames(object@meta)))
})

#' @rdname vars-helpers
#' @export
setMethod("meta_vars", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_vars_fmt(.fmt)
  vars <- lapply(as.list(object), meta_vars)
  if (.fmt == "list") return(vars)
  unique(unlist(vars, use.names = FALSE))
})


#------------------------------------------------------------------------------#
# taxa_vars
#------------------------------------------------------------------------------#

#' @rdname vars-helpers
#' @export
setGeneric("taxa_vars", function(object, .fmt = "list") standardGeneric("taxa_vars"))

#' @rdname vars-helpers
#' @export
setMethod("taxa_vars", "mgnet", function(object, .fmt = "list") {
  if (!has_taxa(object)) return(character(0))
  out <- "taxa_id"
  # add taxa column names only if the taxa slot is present
  if (has_slot(object, "taxa")) out <- c(out, colnames(object@taxa))
  # add comm_id if communities are present
  if (has_comm(object))         out <- c(out, "comm_id")
  unique(out)
})

#' @rdname vars-helpers
#' @export
setMethod("taxa_vars", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_vars_fmt(.fmt)
  vars <- lapply(as.list(object), taxa_vars)
  if (.fmt == "list") return(vars)
  unique(unlist(vars, use.names = FALSE))
})


#------------------------------------------------------------------------------#
# meta_taxa_vars
#------------------------------------------------------------------------------#

#' @rdname vars-helpers
#' @export
setGeneric("meta_taxa_vars", function(object, .fmt = "list") standardGeneric("meta_taxa_vars"))

#' @rdname vars-helpers
#' @export
setMethod("meta_taxa_vars", "mgnet", function(object, .fmt = "list") {
  unique(c(meta_vars(object), taxa_vars(object)))
})

#' @rdname vars-helpers
#' @export
setMethod("meta_taxa_vars", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_vars_fmt(.fmt)
  vars <- lapply(as.list(object), meta_taxa_vars)
  if (.fmt == "list") return(vars)
  unique(c("mgnet", unlist(vars, use.names = FALSE)))
})


#------------------------------------------------------------------------------#
# link_vars
#------------------------------------------------------------------------------#

#' @importFrom igraph edge_attr_names
#' @rdname vars-helpers
#' @export
setGeneric("link_vars", function(object, .fmt = "list") standardGeneric("link_vars"))

#' @rdname vars-helpers
#' @export
setMethod("link_vars", "mgnet", function(object, .fmt = "list") {
  if (!has_netw(object)) return(character(0))
  ea <- igraph::edge_attr_names(object@netw)
  if (length(ea) == 0L)   return("link_id")
  unique(c("link_id", ea))
})

#' @rdname vars-helpers
#' @export
setMethod("link_vars", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_vars_fmt(.fmt)
  vars <- lapply(as.list(object), link_vars)
  if (.fmt == "list") return(vars)
  unique(unlist(vars, use.names = FALSE))
})


#------------------------------------------------------------------------------#
# mgnet_vars
#------------------------------------------------------------------------------#

#' @rdname vars-helpers
#' @export
setGeneric("mgnet_vars", function(object, .fmt = "list") standardGeneric("mgnet_vars"))

#' @rdname vars-helpers
#' @export
setMethod("mgnet_vars", "mgnet", function(object, .fmt = "list") {
  unique(c(meta_vars(object), taxa_vars(object), link_vars(object)))
})

#' @rdname vars-helpers
#' @export
setMethod("mgnet_vars", "mgnets", function(object, .fmt = "list") {
  .fmt <- .match_vars_fmt(.fmt)
  vars <- lapply(as.list(object), mgnet_vars)
  if (.fmt == "list") return(vars)
  unique(c("mgnet", unlist(vars, use.names = FALSE)))
})