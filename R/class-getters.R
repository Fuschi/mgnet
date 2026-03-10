#' @include class-mgnet.R class-mgnets.R class-base-methods.R class-links.R
NULL

#------------------------------------------------------------------------------#
# ABUNDANCE GETTERS
#------------------------------------------------------------------------------#

#' Get abundances from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve abundance data stored in the `abun` slot of an `mgnet` object or
#' across all `mgnet` objects contained in an `mgnets` object.
#'
#' Output can be returned as a matrix, data.frame, or tibble. Optionally,
#' abundances can be aggregated by a taxa-level variable.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of `"matrix"`, `"df"`, or `"tbl"`.
#' @param .var Optional character string giving a taxa-level variable used to
#'   aggregate abundances.
#' @param .fun Aggregation function applied when `.var` is provided.
#' @param .collapse Logical. For `mgnets`, whether to collapse the output into a
#'   single table.
#'
#' @return
#' For `mgnet`, the abundance matrix or its converted representation.
#'
#' For `mgnets`, either:
#' \itemize{
#'   \item a named list of abundance objects, one per `mgnet`, or
#'   \item a single collapsed table if `.collapse = TRUE`.
#' }
#'
#' @details
#' When `.var` is provided, abundances are joined with taxa-level information,
#' aggregated by `sample_id` and `.var`, and then reshaped back to wide format.
#'
#' @export
setGeneric("abun", function(object,
                            .fmt = c("matrix", "df", "tbl"),
                            .var = NULL,
                            .fun = sum,
                            .collapse = FALSE) {
  standardGeneric("abun")
})

#' @rdname abun
setMethod("abun", "mgnet", function(object,
                                    .fmt = c("matrix", "df", "tbl"),
                                    .var = NULL,
                                    .fun = sum,
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (is.null(.var)) {
    if (.fmt == "matrix") return(object@abun)
    if (.fmt == "df")     return(as.data.frame(object@abun))
    if (.fmt == "tbl")    return(tibble::as_tibble(object@abun, rownames = "taxa_id"))
  }
  
  # Ensure .var is a valid column name in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa metadata. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  out <- object@abun |>
    tibble::as_tibble(rownames = "taxa_id") |>
    tidyr::pivot_longer(-"taxa_id", names_to = "sample_id", values_to = "abun") |>
    dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id") |>
    dplyr::group_by(.data$sample_id, .data[[.var]]) |>
    dplyr::summarise(abun = .fun(.data$abun), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "abun")
  
  if (.fmt == "tbl") return(out)
  if (.fmt == "df")  return(as.data.frame(out))
  if (.fmt == "matrix") {
    out <- out |>
      tibble::column_to_rownames(.var) |>
      as.matrix()
    return(out)
  }
})

#' @rdname abun
setMethod("abun", "mgnets", function(object,
                                     .fmt = c("matrix", "df", "tbl"),
                                     .var = NULL,
                                     .fun = sum,
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (.collapse) {
    if (.fmt == "matrix") stop(".fmt = 'matrix' not available when .collapse = TRUE")
    
    out <- object@mgnets |>
      purrr::map(\(x) abun(x, .fmt = "tbl", .var = .var, .fun = .fun)) |>
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
      purrr::list_rbind()
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  purrr::map(object@mgnets, \(x) abun(x, .fmt = .fmt, .var = .var, .fun = .fun))
})


#------------------------------------------------------------------------------#
# RELATIVE ABUNDANCE GETTERS
#------------------------------------------------------------------------------#

#' Get relative abundances from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve relative abundance data stored in the `rela` slot of an `mgnet`
#' object or across all `mgnet` objects contained in an `mgnets` object.
#'
#' Output can be returned as a matrix, data.frame, or tibble. Optionally,
#' relative abundances can be aggregated by a taxa-level variable.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of `"matrix"`, `"df"`, or `"tbl"`.
#' @param .var Optional character string giving a taxa-level variable used to
#'   aggregate relative abundances.
#' @param .fun Aggregation function applied when `.var` is provided.
#' @param .collapse Logical. For `mgnets`, whether to collapse the output into a
#'   single table.
#'
#' @return
#' For `mgnet`, the relative abundance matrix or its converted representation.
#'
#' For `mgnets`, either:
#' \itemize{
#'   \item a named list of relative abundance objects, one per `mgnet`, or
#'   \item a single collapsed table if `.collapse = TRUE`.
#' }
#'
#' @export
setGeneric("rela", function(object,
                            .fmt = c("matrix", "df", "tbl"),
                            .var = NULL,
                            .fun = sum,
                            .collapse = FALSE) {
  standardGeneric("rela")
})

#' @rdname rela
setMethod("rela", "mgnet", function(object,
                                    .fmt = c("matrix", "df", "tbl"),
                                    .var = NULL,
                                    .fun = sum,
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (is.null(.var)) {
    if (.fmt == "matrix") return(object@rela)
    if (.fmt == "df")     return(as.data.frame(object@rela))
    if (.fmt == "tbl")    return(tibble::as_tibble(object@rela, rownames = "taxa_id"))
  }
  
  # Ensure .var is a valid column name in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa metadata. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  out <- object@rela |>
    tibble::as_tibble(rownames = "taxa_id") |>
    tidyr::pivot_longer(-"taxa_id", names_to = "sample_id", values_to = "rela") |>
    dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id") |>
    dplyr::group_by(.data$sample_id, .data[[.var]]) |>
    dplyr::summarise(rela = .fun(.data$rela), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "rela")
  
  if (.fmt == "tbl") return(out)
  if (.fmt == "df")  return(as.data.frame(out))
  if (.fmt == "matrix") {
    out <- out |>
      tibble::column_to_rownames(.var) |>
      as.matrix()
    return(out)
  }
})

#' @rdname rela
setMethod("rela", "mgnets", function(object,
                                     .fmt = c("matrix", "df", "tbl"),
                                     .var = NULL,
                                     .fun = sum,
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (.collapse) {
    if (.fmt == "matrix") stop(".fmt = 'matrix' not available when .collapse = TRUE")
    
    out <- object@mgnets |>
      purrr::map(\(x) rela(x, .fmt = "tbl", .var = .var, .fun = .fun)) |>
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
      purrr::list_rbind()
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  purrr::map(object@mgnets, \(x) rela(x, .fmt = .fmt, .var = .var, .fun = .fun))
})


#------------------------------------------------------------------------------#
# NORMALIZED ABUNDANCE GETTERS
#------------------------------------------------------------------------------#

#' Get normalized abundances from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve normalized abundance data stored in the `norm` slot of an `mgnet`
#' object or across all `mgnet` objects contained in an `mgnets` object.
#'
#' Output can be returned as a matrix, data.frame, or tibble. Optionally,
#' normalized abundances can be aggregated by a taxa-level variable.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of `"matrix"`, `"df"`, or `"tbl"`.
#' @param .var Optional character string giving a taxa-level variable used to
#'   aggregate normalized abundances.
#' @param .fun Aggregation function applied when `.var` is provided.
#' @param .collapse Logical. For `mgnets`, whether to collapse the output into a
#'   single table.
#'
#' @return
#' For `mgnet`, the normalized abundance matrix or its converted representation.
#'
#' For `mgnets`, either:
#' \itemize{
#'   \item a named list of normalized abundance objects, one per `mgnet`, or
#'   \item a single collapsed table if `.collapse = TRUE`.
#' }
#'
#' @export
setGeneric("norm", function(object,
                            .fmt = c("matrix", "df", "tbl"),
                            .var = NULL,
                            .fun = sum,
                            .collapse = FALSE) {
  standardGeneric("norm")
})

#' @rdname norm
setMethod("norm", "mgnet", function(object,
                                    .fmt = c("matrix", "df", "tbl"),
                                    .var = NULL,
                                    .fun = sum,
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (is.null(.var)) {
    if (.fmt == "matrix") return(object@norm)
    if (.fmt == "df")     return(as.data.frame(object@norm))
    if (.fmt == "tbl")    return(tibble::as_tibble(object@norm, rownames = "taxa_id"))
  }
  
  # Ensure .var is a valid column name in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa metadata. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  out <- object@norm |>
    tibble::as_tibble(rownames = "taxa_id") |>
    tidyr::pivot_longer(-"taxa_id", names_to = "sample_id", values_to = "norm") |>
    dplyr::left_join(taxa(object, .fmt = "tbl"), by = "taxa_id") |>
    dplyr::group_by(.data$sample_id, .data[[.var]]) |>
    dplyr::summarise(norm = .fun(.data$norm), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "norm")
  
  if (.fmt == "tbl") return(out)
  if (.fmt == "df")  return(as.data.frame(out))
  if (.fmt == "matrix") {
    out <- out |>
      tibble::column_to_rownames(.var) |>
      as.matrix()
    return(out)
  }
})

#' @rdname norm
setMethod("norm", "mgnets", function(object,
                                     .fmt = c("matrix", "df", "tbl"),
                                     .var = NULL,
                                     .fun = sum,
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  
  if (.collapse) {
    if (.fmt == "matrix") stop(".fmt = 'matrix' not available when .collapse = TRUE")
    
    out <- object@mgnets |>
      purrr::map(\(x) norm(x, .fmt = "tbl", .var = .var, .fun = .fun)) |>
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
      purrr::list_rbind()
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  purrr::map(object@mgnets, \(x) norm(x, .fmt = .fmt, .var = .var, .fun = .fun))
})


#------------------------------------------------------------------------------#
# META GETTERS
#------------------------------------------------------------------------------#

#' Get sample metadata from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve sample metadata stored in the `meta` slot.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of `"df"` or `"tbl"`.
#' @param .empty What to return when no sample metadata are available. One of
#'   `"id"` or `"NULL"`.
#' @param .collapse Logical. For `mgnets`, whether to collapse the output into a
#'   single table.
#'
#' @return
#' For `mgnet`, a data.frame/tibble with sample metadata.
#'
#' For `mgnets`, either a named list of metadata tables or a single collapsed
#' table if `.collapse = TRUE`.
#'
#' @export
setGeneric("meta", function(object,
                            .fmt = c("df", "tbl"),
                            .empty = c("id", "NULL"),
                            .collapse = FALSE) {
  standardGeneric("meta")
})

#' @rdname meta
setMethod("meta", "mgnet", function(object,
                                    .fmt = c("df", "tbl"),
                                    .empty = c("id", "NULL"),
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  .empty <- match.arg(.empty)
  
  if (miss_meta(object)) {
    if (.empty == "NULL") return(NULL)
    
    out <- tibble::tibble(sample_id = sample_id(object))
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  out <- object@meta |>
    tibble::rownames_to_column("sample_id")
  
  if (.fmt == "df") return(as.data.frame(out))
  out
})

#' @rdname meta
setMethod("meta", "mgnets", function(object,
                                     .fmt = c("df", "tbl"),
                                     .empty = c("id", "NULL"),
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  .empty <- match.arg(.empty)
  
  if (length(object) == 0L) {
    if (.fmt == "df") return(data.frame())
    return(tibble::tibble())
  }
  
  if (.collapse) {
    out <- object@mgnets |>
      purrr::map(\(x) meta(x, .fmt = "tbl", .empty = .empty)) |>
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
      purrr::list_rbind()
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  purrr::map(object@mgnets, \(x) meta(x, .fmt = .fmt, .empty = .empty))
})


#------------------------------------------------------------------------------#
# TAXA GETTERS
#------------------------------------------------------------------------------#

#' Get taxa-level information from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve taxa-level information stored in the `taxa` slot and, when present,
#' community assignments stored in the `comm` slot.
#'
#' For a single `mgnet`, the result always includes `taxa_id`. When communities
#' are available, `comm_id` is also included.
#'
#' For `mgnets`, collapsed output always includes both `mgnet` and `taxa_id`.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of `"df"` or `"tbl"`.
#' @param .empty What to return when no taxa-level information is available. One
#'   of `"id"` or `"NULL"`.
#' @param .collapse Logical. For `mgnets`, whether to collapse the output into a
#'   single table.
#'
#' @return
#' For `mgnet`, a data.frame/tibble with taxa-level information.
#'
#' For `mgnets`, either a named list of taxa-level tables or a single collapsed
#' table if `.collapse = TRUE`.
#'
#' @export
setGeneric("taxa", function(object,
                            .fmt = c("df", "tbl"),
                            .empty = c("id", "NULL"),
                            .collapse = FALSE) {
  standardGeneric("taxa")
})

#' @rdname taxa
setMethod("taxa", "mgnet", function(object,
                                    .fmt = c("df", "tbl"),
                                    .empty = c("id", "NULL"),
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  .empty <- match.arg(.empty)
  
  # No taxa metadata and no community assignments
  if (miss_taxa(object) && miss_comm(object)) {
    if (.empty == "NULL") return(NULL)
    
    out <- tibble::tibble(taxa_id = taxa_id(object))
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  # Taxa metadata only
  if (!miss_taxa(object) && miss_comm(object)) {
    out <- object@taxa |>
      tibble::rownames_to_column("taxa_id")
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  # Community assignments only
  if (miss_taxa(object) && !miss_comm(object)) {
    out <- comm_id(object, .fmt = "tbl")
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  # Both taxa metadata and community assignments
  out <- object@taxa |>
    tibble::rownames_to_column("taxa_id") |>
    dplyr::left_join(comm_id(object, .fmt = "tbl"), by = "taxa_id")
  
  if (.fmt == "df") return(as.data.frame(out))
  out
})

#' @rdname taxa
setMethod("taxa", "mgnets", function(object,
                                     .fmt = c("df", "tbl"),
                                     .empty = c("id", "NULL"),
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt)
  .empty <- match.arg(.empty)
  
  if (length(object) == 0L) {
    if (.fmt == "df") return(data.frame())
    return(tibble::tibble())
  }
  
  if (.collapse) {
    out <- object@mgnets |>
      purrr::map(\(x) taxa(x, .fmt = "tbl", .empty = .empty)) |>
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
      purrr::list_rbind()
    
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  purrr::map(object@mgnets, \(x) taxa(x, .fmt = .fmt, .empty = .empty))
})


#------------------------------------------------------------------------------#
# NETWORK GETTERS
#------------------------------------------------------------------------------#

#' Get network objects from `mgnet` and `mgnets`
#'
#' @description
#' Retrieve the network stored in the `netw` slot.
#'
#' If taxa-level information is available, it is added as vertex attributes
#' before returning the graph.
#'
#' When `selected = TRUE`, only currently selected links are retained. If the
#' current selection is empty, the returned graph contains zero edges while
#' preserving all vertices.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param selected Logical. Whether to keep only selected links.
#'
#' @return
#' For `mgnet`, an `igraph` object.
#'
#' For `mgnets`, a named list of `igraph` objects.
#'
#' @export
setGeneric("netw", function(object, selected = FALSE) {
  standardGeneric("netw")
})

#' @rdname netw
setMethod("netw", "mgnet", function(object, selected = FALSE) {
  g <- object@netw
  
  if (!miss_metataxa(object)) {
    tx <- taxa(object, .fmt = "tbl", .empty = "id")
    
    for (nm in names(tx)[-1]) {
      vals <- tx[[nm]]
      names(vals) <- tx$taxa_id
      g <- igraph::set_vertex_attr(g, nm, value = vals[igraph::V(g)$name])
    }
  }
  
  if (selected && are_selected_links(object)) {
    keep <- igraph::edge_attr(g, "link_id") %in% get_selected_links(object)
    g <- igraph::subgraph_from_edges(g, eids = which(keep), delete.vertices = FALSE)
  }
  
  g
})

#' @rdname netw
setMethod("netw", "mgnets", function(object, selected = FALSE) {
  purrr::map(object@mgnets, \(x) netw(x, selected = selected))
})


#------------------------------------------------------------------------------#
# LINK GETTERS
#------------------------------------------------------------------------------#

#' Get link-level data from `mgnet` and `mgnets` objects
#'
#' @description
#' Retrieve edge tables from the network stored in an `mgnet` or `mgnets`
#' object.
#'
#' By default, only selected links are returned.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param selected Logical. Whether to return only selected links. Default is
#'   `TRUE`.
#' @param .fmt For `mgnets`, output format. One of `"list"` or `"tbl"`.
#'
#' @return
#' For `mgnet`, a data.frame of links.
#'
#' For `mgnets`, either a named list of link tables or a single collapsed tibble
#' with an additional `mgnet` column.
#'
#' @details
#' This method calls [netw()] and then converts the graph(s) to edge tables with
#' `igraph::as_data_frame(..., what = "edges")`.
#'
#' @export
setGeneric("link", function(object, selected = TRUE, .fmt = c("list", "tbl")) {
  standardGeneric("link")
})

#' @rdname link
setMethod("link", "mgnet", function(object, selected = TRUE, .fmt = c("list", "tbl")) {
  if (miss_netw(object)) {
    cli::cli_abort("No network available for {.cls mgnet} object.")
  }
  
  igraph::as_data_frame(netw(object, selected = selected), what = "edges")
})

#' @rdname link
setMethod("link", "mgnets", function(object, selected = TRUE, .fmt = c("list", "tbl")) {
  .fmt <- match.arg(.fmt)
  
  if (miss_netw(object, .fmt = "any")) {
    cli::cli_abort("No network available in at least one element of the {.cls mgnets} object.")
  }
  
  out <- purrr::map(object@mgnets, \(x) link(x, selected = selected))
  
  if (.fmt == "list") return(out)
  
  out |>
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
    purrr::list_rbind()
})