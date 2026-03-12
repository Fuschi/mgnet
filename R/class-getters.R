#' @include class-mgnet.R class-mgnets.R class-base-methods.R class-links.R
NULL

# Internal helpers shared across getters
#------------------------------------------------------------------------------#
.as_abundance_output <- function(x, .fmt) {
  switch(.fmt,
         mat = x,
         df  = as.data.frame(x),
         tbl = tibble::as_tibble(x, rownames = "sample_id"))
}

.aggregate_abundance_by_var <- function(object, x, value_name, .var, .fun) {
  data_frame <- x %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = value_name)
  
  taxa_info <- taxa(object, .fmt = "tbl") %>%
    dplyr::select(taxa_id, !!rlang::sym(.var))
  
  data_frame %>%
    dplyr::left_join(taxa_info, by = "taxa_id") %>%
    dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
    dplyr::summarise(!!rlang::sym(value_name) := .fun(.data[[value_name]]), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = !!rlang::sym(value_name)) %>%
    dplyr::arrange(match(sample_id, sample_id(object))) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
}

.get_abundance_slot <- function(object, slot_name, .fmt = "mat", .var = NULL, .fun = sum) {
  if (!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if (!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  x <- methods::slot(object, slot_name)
  
  if (length(x) == 0) {
    return(
      switch(.fmt,
             mat = matrix(nrow = 0, ncol = 0),
             df  = data.frame(),
             tbl = tibble::tibble())
    )
  }
  
  if (is.null(.var)) {
    return(.as_abundance_output(x, .fmt))
  }
  
  available <- taxa_vars(object)
  if (!(.var %in% available)) {
    stop(
      ".var must be present in the taxa-level information. Available choices are: ",
      paste(available, collapse = ", ")
    )
  }
  
  result_matrix <- .aggregate_abundance_by_var(object, x, slot_name, .var, .fun)
  .as_abundance_output(result_matrix, .fmt)
}

.map_mgnets <- function(object, fun) {
  purrr::map(object@mgnets, fun)
}

.collapse_mgnets_tbl <- function(object, fun) {
  object@mgnets |>
    purrr::map(fun) |>
    purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) |>
    purrr::list_rbind()
}

.empty_collection_output <- function(.fmt) {
  if (.fmt == "df") return(data.frame())
  tibble::tibble()
}

.format_meta_output <- function(object, .fmt) {
  switch(.fmt,
         df  = object@meta,
         tbl = tibble::as_tibble(object@meta, rownames = "sample_id"))
}

.format_taxa_output <- function(object, .fmt) {
  if (miss_taxa(object) && miss_comm(object)) {
    out <- tibble::tibble(taxa_id = taxa_id(object))
    return(if (.fmt == "df") as.data.frame(out) else out)
  }
  
  if (!miss_taxa(object) && miss_comm(object)) {
    return(
      switch(.fmt,
             df  = object@taxa,
             tbl = tibble::as_tibble(object@taxa, rownames = "taxa_id"))
    )
  }
  
  if (miss_taxa(object) && !miss_comm(object)) {
    out <- comm_id(object, .fmt = "tbl")
    return(if (.fmt == "df") as.data.frame(out) else out)
  }
  
  out <- tibble::as_tibble(object@taxa, rownames = "taxa_id") |>
    dplyr::left_join(comm_id(object, .fmt = "tbl"), by = "taxa_id")
  
  switch(.fmt,
         df  = tibble::column_to_rownames(out, "taxa_id"),
         tbl = out)
}

#------------------------------------------------------------------------------#
# ABUNDANCE / RELATIVE ABUNDANCE / NORMALIZED ABUNDANCE
#------------------------------------------------------------------------------#

#' Retrieve abundance-like matrices from `mgnet` and `mgnets`
#'
#' @description
#' Retrieve and optionally aggregate abundance-like matrices stored in an
#' `mgnet` or `mgnets` object.
#'
#' These methods support:
#' \itemize{
#'   \item raw abundance via [abun()],
#'   \item relative abundance via [rela()],
#'   \item normalized abundance via [norm()].
#' }
#'
#' If `.var` is provided, values are aggregated according to a taxa-level
#' variable returned by [taxa_vars()]. Aggregation is performed independently
#' within each sample.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt Output format. One of:
#' \itemize{
#'   \item `"mat"`: return a matrix
#'   \item `"df"`: return a data.frame
#'   \item `"tbl"`: return a tibble
#' }
#' @param .var Optional character string giving a taxa-level variable used for
#'   grouping and aggregation. If omitted, the original matrix is returned.
#' @param .fun A function used to aggregate values within each sample and group.
#'   Defaults to `sum`.
#'
#' @return
#' For an `mgnet` object:
#' \itemize{
#'   \item if `.var` is `NULL`, the original matrix in the requested format;
#'   \item otherwise, an aggregated matrix/table with rows corresponding to
#'   `sample_id` and columns corresponding to the levels of `.var`.
#' }
#'
#' For an `mgnets` object:
#' \itemize{
#'   \item a named list containing the corresponding result for each contained
#'   `mgnet`.
#' }
#'
#' @details
#' The abundance-like matrices are assumed to have:
#' \itemize{
#'   \item rows = `sample_id`
#'   \item columns = `taxa_id`
#' }
#'
#' When `.var` is provided, the matrix is converted to long format, joined with
#' taxa-level information, aggregated by `sample_id` and `.var`, and then
#' reshaped back to wide format.
#'
#' @name abundance-getters
NULL


#------------------------------------------------------------------------------#
# ABUNDANCE
#------------------------------------------------------------------------------#

#' @rdname abundance-getters
#' @export
#' @name abun
#' @aliases abun,mgnet-method abun,mgnets-method
setGeneric("abun", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  standardGeneric("abun")
})

#' @rdname abundance-getters
#' @export
setMethod("abun", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .get_abundance_slot(object, "abun", .fmt = .fmt, .var = .var, .fun = .fun)
})

#' @rdname abundance-getters
#' @export
setMethod("abun", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  .map_mgnets(object, \(x) abun(object = x, .var = .var, .fmt = .fmt, .fun = .fun))
})


#------------------------------------------------------------------------------#
# RELATIVE ABUNDANCE
#------------------------------------------------------------------------------#

#' @rdname abundance-getters
#' @export
#' @name rela
#' @aliases rela,mgnet-method rela,mgnets-method
setGeneric("rela", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  standardGeneric("rela")
})

#' @rdname abundance-getters
#' @export
setMethod("rela", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .get_abundance_slot(object, "rela", .fmt = .fmt, .var = .var, .fun = .fun)
})

#' @rdname abundance-getters
#' @export
setMethod("rela", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  .map_mgnets(object, \(x) rela(object = x, .var = .var, .fmt = .fmt, .fun = .fun))
})


#------------------------------------------------------------------------------#
# NORMALIZED ABUNDANCE
#------------------------------------------------------------------------------#

#' @rdname abundance-getters
#' @export
#' @name norm
#' @aliases norm,mgnet-method norm,mgnets-method
setGeneric("norm", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  standardGeneric("norm")
})

#' @rdname abundance-getters
#' @export
setMethod("norm", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .get_abundance_slot(object, "norm", .fmt = .fmt, .var = .var, .fun = .fun)
})

#' @rdname abundance-getters
#' @export
setMethod("norm", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  .map_mgnets(object, \(x) norm(object = x, .var = .var, .fmt = .fmt, .fun = .fun))
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
                            .fmt = "df",
                            .empty = "id",
                            .collapse = FALSE) {
  standardGeneric("meta")
})

#' @rdname meta
setMethod("meta", "mgnet", function(object,
                                    .fmt = "df",
                                    .empty = "id",
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  .empty <- match.arg(.empty, c("id", NULL))
  
  if (miss_meta(object)) {
    if (.empty == "NULL") return(NULL)
    
    out <- tibble::tibble(sample_id = sample_id(object))
    if (.fmt == "df") return(as.data.frame(out))
    return(out)
  }
  
  .format_meta_output(object, .fmt)
})

#' @rdname meta
setMethod("meta", "mgnets", function(object,
                                     .fmt = "df",
                                     .empty = "id",
                                     .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  .empty <- match.arg(.empty, c("id", NULL))
  
  if (length(object) == 0L) {
    return(.empty_collection_output(.fmt))
  }
  
  if (.collapse) {
    return(.collapse_mgnets_tbl(object, \(x) meta(x, .fmt = "tbl", .empty = .empty)))
  }
  
  .map_mgnets(object, \(x) meta(x, .fmt = .fmt, .empty = .empty))
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
                            .fmt = "df",
                            .empty = "id",
                            .collapse = FALSE) {
  standardGeneric("taxa")
})

#' @rdname taxa
setMethod("taxa", "mgnet", function(object,
                                    .fmt = "df",
                                    .empty = "id",
                                    .collapse = FALSE) {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  .empty <- match.arg(.empty, c("id", NULL))
  
  if (miss_taxa(object) && miss_comm(object) && .empty == "NULL") return(NULL)
  .format_taxa_output(object, .fmt)
})

#' @rdname taxa
setMethod("taxa", "mgnets", function(object,
                                     .fmt = "df",
                                     .empty = "id",
                                     .collapse = FALSE) {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  .empty <- match.arg(.empty, c("id", NULL))
  
  if (length(object) == 0L) {
    return(.empty_collection_output(.fmt))
  }
  
  if (.collapse) {
    return(.collapse_mgnets_tbl(object, \(x) taxa(x, .fmt = "tbl", .empty = .empty)))
  }
  
  .map_mgnets(object, \(x) taxa(x, .fmt = .fmt, .empty = .empty))
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
#' @param selected Logical (default TRUE). Whether to keep only selected links.
#'
#' @return
#' For `mgnet`, an `igraph` object.
#'
#' For `mgnets`, a named list of `igraph` objects.
#'
#' @export
setGeneric("netw", function(object, selected = TRUE) {
  standardGeneric("netw")
})

#' @rdname netw
setMethod("netw", "mgnet", function(object, selected = TRUE) {
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
setMethod("netw", "mgnets", function(object, selected = TRUE) {
  .map_mgnets(object, \(x) netw(x, selected = selected))
})

#------------------------------------------------------------------------------#
# COMM GETTERS
#------------------------------------------------------------------------------#

#' Access community detection results
#'
#' Retrieve the community detection object stored in the `comm` slot of an
#' `mgnet` object, or from each element of an `mgnets` object.
#'
#' @param object An `mgnet` or `mgnets` object.
#'
#' @return
#' For `mgnet`, the community detection object stored in the `comm` slot.
#'  
#' For `mgnets`, a named list of community detection objects.
#'
#' @export
#' @name comm
#' @aliases comm,mgnet-method comm,mgnets-method
setGeneric("comm", function(object) standardGeneric("comm"))

setMethod("comm", "mgnet", function(object) {
  object@comm
})

setMethod("comm", "mgnets", function(object) {
  .map_mgnets(object, \(x) x@comm)
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
  
  out <- .map_mgnets(object, \(x) link(x, selected = selected))
  
  if (.fmt == "list") return(out)
  
  .collapse_mgnets_tbl(object, \(x) link(x, selected = selected))
})
