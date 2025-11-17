#------------------------------------------------------------------------------#
#' Join and Convert Abundance Data to Long Format
#'
#' This internal function creates a grid of `sample_id` and `taxa_id` combinations 
#' and joins all required abundance variables in long format. It converts abundance data 
#' (e.g., `abun`, `rela`, `norm`) from wide to long format for easier manipulation.
#'
#' @param object An `mgnet` object.
#' @param needed_keys A list containing the abundance keys to be used (e.g., `abun`, `rela`, `norm`).
#' 
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom dplyr left_join left_join
#' @importFrom methods slot
#' @importFrom tibble as_tibble
#' @importFrom purrr map imap list_rbind
#'
#' @return A tibble with `sample_id`, `taxa_id`, and the abundance variables in long format.
#' @keywords internal
long_abundance_join <- function(object, needed_keys) {
  
  if(inherits(object, "mgnet")){
    
    if (miss_sample(object) || miss_taxa(object)) {
      cli::cli_abort("No sample or taxa available in the {.cls mgnet} object.")}
    if(length(needed_keys) == 0) return(tibble::tibble())
    
    long_abun <- tidyr::expand_grid(
      sample_id = sample_id(object),
      taxa_id = taxa_id(object))
    
    for (abundance_key in needed_keys) {
      long_abun <- long_abun %>%
        dplyr::left_join(
          methods::slot(object, abundance_key) %>%
            tibble::as_tibble(rownames = "sample_id") %>%
            tidyr::pivot_longer(-sample_id,
                                names_to = "taxa_id",
                                values_to = abundance_key),
          by = c("sample_id", "taxa_id")
        )}
    
    return(long_abun)
    
  } else if (inherits(object, "mgnets")){
    
    if(miss_sample(object, "any") || miss_taxa(object, "any")) stop("Error: No sample or taxa available in at least one of the mgnet objects.")
    
    long_abun <- purrr::map(object, long_abundance_join, needed_keys = needed_keys) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
    
    return(long_abun)
    
  } else {
    
    cli::cli_abort(c(
      "x" = "Unsupported object class.",
      "i" = "Expected {.cls mgnet} or {.cls mgnets}, got: {.cls {class(object)}}."
    ))
    
  }
}