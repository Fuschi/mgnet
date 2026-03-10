#' @include class-mgnet.R class-mgnets.R class-base-methods.R class-links.R
NULL

# ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Abundance Data 
#'
#' @description
#' Retrieves and aggregates abundance data for an `mgnet` or `mgnets` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnets` object containing abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original abundance data without any aggregation.
#' @param .fun A function to specify how the abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnets` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name abun
#' @aliases abun,mgnet-method abun,mgnets-method
setGeneric("abun", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("abun"))

setMethod("abun", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty abundance data
  if(length(object@abun) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@abun,
           df = as.data.frame(object@abun),
           tbl = tibble::as_tibble(object@abun, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% taxa_vars(object))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@abun %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "abun")
    
    taxonmgnet_info <- taxa(object, .fmt = "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonmgnet_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(abun = .fun(abun), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "abun") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("abun", "mgnets", function(object,  .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@mgnets, function(x) abun(object = x, .var = .var, 
                                                   .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# RELATIVE ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Relative Abundance Data Based on Taxonmgnet Information
#'
#' @description
#' Retrieves and aggregates relative abundance data for an `mgnet` or `mgnets` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original relative abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain relative abundances at different taxonmgnet ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnets` object containing relative abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the relative abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original relative abundance data without any aggregation.
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns relative abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnets` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name rela
#' @aliases rela,mgnet-method rela,mgnets-method
setGeneric("rela", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("rela"))

setMethod("rela", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty rela data
  if(length(object@rela) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@rela,
           df = as.data.frame(object@rela),
           tbl = tibble::as_tibble(object@rela, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% colnames(object@taxa))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@rela %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "rela")
    
    taxonmgnet_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonmgnet_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(rela = .fun(rela), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "rela") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("rela", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@mgnets, function(x) rela(object = x, .var = .var, 
                                                   .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# NORMALIZED ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Normalized Abundance Data Based on Taxonmgnet Information
#'
#' @description
#' Retrieves and aggregates normalized abundance data for an `mgnet` or `mgnets` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original normalized abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain normalized abundances at different taxonmgnet ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnets` object containing normalized abundance and informational data.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the normalized abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original normalized abundance data without any aggregation.
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns normalized abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnets` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name norm
#' @aliases norm,mgnet-method norm,mgnets-method
setGeneric("norm", function(object, .fmt = "mat", .var = NULL, .fun = sum) standardGeneric("norm"))

setMethod("norm", "mgnet", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  # Checks
  if(!is.null(.var) && !is.character(.var)) {
    stop(".var must be a character string specifying the column name.")
  }
  
  if(!is.function(.fun)) {
    stop(".fun must be a function.")
  }
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  # Handling empty norm data
  if(length(object@norm) == 0) {
    switch(.fmt,
           mat = matrix(nrow = 0, ncol = 0),
           df = data.frame(),
           tbl = tibble::tibble())
  } else if(is.null(.var)) {
    # Return data in specified format without aggregation
    switch(.fmt,
           mat = object@norm,
           df = as.data.frame(object@norm),
           tbl = tibble::as_tibble(object@norm, rownames = "sample_id"))
  } else {
    # Ensure .var is a valid column name in taxa
    if(!(.var %in% colnames(object@taxa))) {
      stop(".var must be present in the taxa columns. Available choices are: ", 
           paste(colnames(object@taxa), collapse = ", "))
    }
    
    # Prepare data for aggregation
    data_frame <- object@norm %>% 
      tibble::as_tibble(rownames = "sample_id") %>%
      tidyr::pivot_longer(cols = -sample_id, names_to = "taxa_id", values_to = "norm")
    
    taxonmgnet_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonmgnet_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(.var)) %>%
      dplyr::summarise(norm = .fun(norm), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(.var), values_from = "norm") %>%
      dplyr::arrange(match(sample_id, sample_id(object)))
    
    result_matrix <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result_matrix)},
           df  = {return(as.data.frame(result_matrix))},
           tbl = {return(tibble::as_tibble(result_matrix, rownames="sample_id"))}
    )
  }
})

setMethod("norm", "mgnets", function(object, .fmt = "mat", .var = NULL, .fun = sum) {
  
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object, function(x) norm(object = x, .var = .var, .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
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