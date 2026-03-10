#' @include class-mgnet.R class-mgnets.R class-base-methods.R class-links.R
NULL

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
  
  # Checks
  if (!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if (!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handle empty abundance slot
  if (length(object@abun) == 0) {
    return(
      switch(.fmt,
             mat = matrix(nrow = 0, ncol = 0),
             df  = data.frame(),
             tbl = tibble::tibble())
    )
  }
  
  # No aggregation requested
  if (is.null(.var)) {
    return(
      switch(.fmt,
             mat = object@abun,
             df  = as.data.frame(object@abun),
             tbl = tibble::as_tibble(object@abun, rownames = "sample_id"))
    )
  }
  
  # Ensure .var is available in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa-level information. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  # Prepare long table: rows = sample_id, columns = taxa_id
  data_frame <- object@abun %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "abun")
  
  taxa_info <- taxa(object, .fmt = "tbl") %>%
    dplyr::select(taxa_id, !!rlang::sym(.var))
  
  # Aggregate within each sample and taxa grouping
  aggregated_data <- data_frame %>%
    dplyr::left_join(taxa_info, by = "taxa_id") %>%
    dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
    dplyr::summarise(abun = .fun(abun), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "abun") %>%
    dplyr::arrange(match(sample_id, sample_id(object)))
  
  result_matrix <- aggregated_data %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  
  switch(.fmt,
         mat = result_matrix,
         df  = as.data.frame(result_matrix),
         tbl = tibble::as_tibble(result_matrix, rownames = "sample_id"))
})

#' @rdname abundance-getters
#' @export
setMethod("abun", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  sapply(
    object@mgnets,
    function(x) abun(object = x, .var = .var, .fmt = .fmt, .fun = .fun),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
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
  
  # Checks
  if (!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if (!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handle empty relative abundance slot
  if (length(object@rela) == 0) {
    return(
      switch(.fmt,
             mat = matrix(nrow = 0, ncol = 0),
             df  = data.frame(),
             tbl = tibble::tibble())
    )
  }
  
  # No aggregation requested
  if (is.null(.var)) {
    return(
      switch(.fmt,
             mat = object@rela,
             df  = as.data.frame(object@rela),
             tbl = tibble::as_tibble(object@rela, rownames = "sample_id"))
    )
  }
  
  # Ensure .var is available in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa-level information. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  # Prepare long table: rows = sample_id, columns = taxa_id
  data_frame <- object@rela %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "rela")
  
  taxa_info <- taxa(object, .fmt = "tbl") %>%
    dplyr::select(taxa_id, !!rlang::sym(.var))
  
  # Aggregate within each sample and taxa grouping
  aggregated_data <- data_frame %>%
    dplyr::left_join(taxa_info, by = "taxa_id") %>%
    dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
    dplyr::summarise(rela = .fun(rela), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "rela") %>%
    dplyr::arrange(match(sample_id, sample_id(object)))
  
  result_matrix <- aggregated_data %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  
  switch(.fmt,
         mat = result_matrix,
         df  = as.data.frame(result_matrix),
         tbl = tibble::as_tibble(result_matrix, rownames = "sample_id"))
})

#' @rdname abundance-getters
#' @export
setMethod("rela", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  sapply(
    object@mgnets,
    function(x) rela(object = x, .var = .var, .fmt = .fmt, .fun = .fun),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
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
  
  # Checks
  if (!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if (!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handle empty normalized abundance slot
  if (length(object@norm) == 0) {
    return(
      switch(.fmt,
             mat = matrix(nrow = 0, ncol = 0),
             df  = data.frame(),
             tbl = tibble::tibble())
    )
  }
  
  # No aggregation requested
  if (is.null(.var)) {
    return(
      switch(.fmt,
             mat = object@norm,
             df  = as.data.frame(object@norm),
             tbl = tibble::as_tibble(object@norm, rownames = "sample_id"))
    )
  }
  
  # Ensure .var is available in taxa-level information
  if (!(.var %in% taxa_vars(object))) {
    stop(
      ".var must be present in the taxa-level information. Available choices are: ",
      paste(taxa_vars(object), collapse = ", ")
    )
  }
  
  # Prepare long table: rows = sample_id, columns = taxa_id
  data_frame <- object@norm %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "norm")
  
  taxa_info <- taxa(object, .fmt = "tbl") %>%
    dplyr::select(taxa_id, !!rlang::sym(.var))
  
  # Aggregate within each sample and taxa grouping
  aggregated_data <- data_frame %>%
    dplyr::left_join(taxa_info, by = "taxa_id") %>%
    dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
    dplyr::summarise(norm = .fun(norm), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "norm") %>%
    dplyr::arrange(match(sample_id, sample_id(object)))
  
  result_matrix <- aggregated_data %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  
  switch(.fmt,
         mat = result_matrix,
         df  = as.data.frame(result_matrix),
         tbl = tibble::as_tibble(result_matrix, rownames = "sample_id"))
})

#' @rdname abundance-getters
#' @export
setMethod("norm", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  sapply(
    object@mgnets,
    function(x) norm(object = x, .var = .var, .fmt = .fmt, .fun = .fun),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
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
  purrr::map(object@mgnets, \(x) netw(x, selected = selected))
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
  lapply(object@mgnets, function(x) x@comm)
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