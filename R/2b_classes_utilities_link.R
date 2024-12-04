#' Select Links in mgnet Object Based on Conditions
#'
#' This method filters links in the `mgnet` object according to specified conditions and stores the indices
#' of the selected links as an attribute of the object.
#'
#' @param object An `mgnet` object.
#' @param ... Conditions to be applied for selecting links, passed as expressions.
#' @param .suffix Character vector specifying suffixes for identifying nodes in links.
#' @return Returns the `mgnet` object with the `selected_links` attribute set to the indices of the filtered links.
#' @importFrom rlang enquos
#' @importFrom dplyr  mutate filter pull
#' @export
setGeneric("select_link", function(object, ..., .suffix = c("_1", "_2")) standardGeneric("select_link"))

setMethod("select_link", "mgnet", function(object, ..., .suffix = c("_1", "_2")) {

  quosures <- rlang::enquos(...)
  object <- deselect_links(object)
  edges_df <- gather_link(object, .suffix = .suffix)
  
  filtered_links <- edges_df %>%
    dplyr::mutate(link_id = row_number()) %>%
    dplyr::filter(!!!quosures) %>%
    dplyr::pull(link_id)
  
  # Store the filtered links as an attribute
  attr(object, "selected_links") <- filtered_links
  return(object)
})

setMethod("select_link", "mgnetList", function(object, ..., .suffix = c("_1", "_2")) {
  
  quosures <- rlang::enquos(...)
  object <- deselect_links(object)
  edges_df <- gather_link(object, .suffix = .suffix)
  
  filtered_links <- edges_df %>%
    dplyr::group_by(mgnet) %>%
    dplyr::mutate(link_id = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!!quosures) %>%
    dplyr::select(mgnet, link_id)
  
  filtered_links_splitted <-  filtered_links %>%
    base::split(.[, "mgnet"]) %>%
    purrr::map(\(x) pull(x, link_id))
  
  # Store the filtered links as an attribute
  for(i in 1:length(object)){
    attr(object[[i]], "selected_links") <- filtered_links_splitted[[i]]
  }
  return(object)
  
  
})

#' Get Selected Links from an mgnet Object
#'
#' Retrieves the indices of the links that were previously selected using `select_link`.
#'
#' @param object An `mgnet` object.
#' @return Numeric vector of indices representing the selected links.
#' @export
setGeneric("get_selected_links", function(object) standardGeneric("get_selected_links"))

setMethod("get_selected_links", "mgnet", function(object) {
  return(attr(object, "selected_links"))
})

setMethod("get_selected_links", "mgnetList", function(object) {
  sapply(object, get_selected_links, USE.NAMES = TRUE, simplify = FALSE)
})

#' Deselect Links in an mgnet Object
#'
#' Clears the `selected_links` attribute from the `mgnet` object, effectively removing the selection of links.
#' This method resets the state of link selection within the object, allowing for new operations without the constraints of the previous link selection.
#'
#' @param object An `mgnet` object.
#' @return The `mgnet` object without the `selected_links` attribute, reflecting no current link selection.
#' @export
setGeneric("deselect_links", function(object) standardGeneric("deselect_links"))

setMethod("deselect_links", "mgnet", function(object) {
  attr(object, "selected_links") <- NULL
  return(object)
})

setMethod("deselect_links", "mgnetList", function(object) {
  for(i in seq_along(object)) object[[i]] <- deselect_links(object[[i]])
  return(object)
})