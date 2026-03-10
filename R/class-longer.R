#' @include class-mgnet.R class-mgnets.R class-base-methods.R
NULL

#' Convert an `mgnet` or `mgnets` object to a long format table
#'
#' @description
#' `mgnet_longer()` converts an `mgnet` object into a long-format tibble
#' containing all `sample_id` × `taxa_id` combinations, together with any
#' available abundance matrices (`abun`, `rela`, `norm`) and associated
#' metadata (e.g. `meta`, `taxa`, `comm_id`).
#'
#' For `mgnets`, the method applies the same transformation to each contained
#' `mgnet` object and combines the results into a single tibble with an
#' additional `mgnet` column.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' A tibble in long format.
#'
#' For `mgnet`, the tibble contains:
#' \itemize{
#'   \item `sample_id`
#'   \item `taxa_id`
#'   \item optional abundance columns (`abun`, `rela`, `norm`)
#'   \item sample metadata from `meta(object)`
#'   \item taxa-level information from `taxa(object)`
#' }
#'
#' For `mgnets`, the tibble additionally contains:
#' \itemize{
#'   \item `mgnet`: the name of the originating `mgnet` object
#' }
#'
#' @details
#' The output includes all sample–taxa combinations, even when abundance values
#' are missing, so abundance columns may contain `NA`.
#'
#' @name mgnet_longer
#' @rdname mgnet_longer
#' @export
setGeneric("mgnet_longer", function(object) standardGeneric("mgnet_longer"))

#' @rdname mgnet_longer
#' @export
setMethod("mgnet_longer", "mgnet", function(object) {
  
  # Basic availability checks
  if (miss_sample(object) || miss_taxa(object)) {
    cli::cli_abort("Error: No sample or taxa available.")
  }
  
  # Start with the full sample × taxa grid
  long_mgnet <- tidyr::expand_grid(
    sample_id = sample_id(object),
    taxa_id   = taxa_id(object)
  )
  
  # Add abundances if available
  if (has_slot(object, "abun")) {
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        abun(object, .fmt = "tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "abun"),
        by = c("sample_id", "taxa_id")
      )
  }
  
  # Add relative abundances if available
  if (has_slot(object, "rela")) {
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        rela(object, .fmt = "tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "rela"),
        by = c("sample_id", "taxa_id")
      )
  }
  
  # Add normalized abundances if available
  if (has_slot(object, "norm")) {
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        norm(object, .fmt = "tbl") %>%
          tidyr::pivot_longer(-sample_id, names_to = "taxa_id", values_to = "norm"),
        by = c("sample_id", "taxa_id")
      )
  }
  
  # Add sample metadata if available
  if (has_slot(object, "meta")) {
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        meta(object, .fmt = "tbl"),
        by = "sample_id"
      )
  }
  
  # Add taxa-level information if available
  if (has_metataxa(object)) {
    long_mgnet <- long_mgnet %>%
      dplyr::left_join(
        taxa(object, .fmt = "tbl"),
        by = "taxa_id"
      )
  }
  
  long_mgnet
})


#' @rdname mgnet_longer
#' @export
setMethod("mgnet_longer", "mgnets", function(object) {
  
  if (miss_sample(object, "any") || miss_taxa(object, "any")) {
    cli::cli_abort("Error: No sample or taxa available in at least one mgnet object.")
  }
  
  purrr::map(object, mgnet_longer) %>%
    purrr::imap(\(x, y) {
      dplyr::mutate(x, mgnet = y, .before = 1)
    }) %>%
    purrr::list_rbind()
})