#'@include class-mgnet.R class-mgnets.R
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


# META
#------------------------------------------------------------------------------#
#' Get Sample Metadata Information
#'
#' Retrieves the sample information stored in the `meta` slot of an `mgnet` object
#' or for each `mgnet` object within an `mgnets`, with the option to format the output as
#' a `data.frame`, `tibble`.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, the row names of
#'          the abundance matrix are converted into a new column named `sample_id`, ensuring alignment
#'          with the reserved keyword in `mgnet-class`.
#' @param .collapse Logical, only for `mgnets`. If TRUE, ignore `.fmt` and return
#'   a single tibble with all metas row-bound together and an `mgnet` column  
#'   indicating the source object. Default: FALSE.
#'          
#' @return The content of the `meta` slot for `mgnet` or a list of such contents for `mgnets`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name meta
#' @aliases meta,mgnet-method meta,mgnets-method
setGeneric("meta", function(object, .fmt = "df", .collapse = FALSE) standardGeneric("meta"))

setMethod("meta", "mgnet", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if( length(object@meta) == 0 ){
    
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
    
  } else {
    
    switch(.fmt,
           df  = {return(object@meta)},
           tbl = {return(tibble::as_tibble(object@meta, rownames="sample_id"))})
    
  }
})

setMethod("meta", "mgnets", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if (.collapse) {
    meta_collapsed <- object %>% 
      purrr::map(\(x) meta(x, .fmt = "tbl")) %>% 
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>% 
      purrr::list_rbind()
      return(meta_collapsed)
  }
  
  if(.fmt == "df") {
    
    if(length(object) == 0) return(data.frame())
    return(purrr::map(object, function(x) meta(x, .fmt = "df")))
    
  } else if(.fmt == "tbl") {
    
    if(length(object) == 0) return(tibble::tibble())
    return(purrr::map(object, function(x) meta(x, .fmt = "tbl")))
    
  }
  
})


# TAXA
#------------------------------------------------------------------------------#
#' Get Taxa Metadata Information
#'
#' Retrieves the taxa information stored in the `taxa` slot of an `mgnet` object
#' or for each `mgnet` object within an `mgnets`, with the option to format the output as
#' a `data.frame`, `tibble`, or a combined `tibble` for multiple mgnet objects.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, the the row names of
#'          the abundance matrix are converted into a new column named `taxa_id`, ensuring alignment
#'          with the reserved keyword in `mgnet-class`.
#' @param .collapse Logical, only for `mgnets`. If TRUE, ignore `.fmt` and return
#'   a single tibble with all taxas row-bound together and an `mgnet` column  
#'   indicating the source object. Default: FALSE.
#'          
#' @return The content of the `taxa` slot for `mgnet` or a list of such contents for `mgnets`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name taxa
#' @aliases taxa,mgnet-method taxa,mgnets-method
setGeneric("taxa", function(object, .fmt = "df", .collapse = FALSE) standardGeneric("taxa"))

setMethod("taxa", "mgnet", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if(length(object@taxa) == 0 && length(object@comm) == 0){
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  }
  
  if(length(object@taxa) != 0 && length(object@comm) == 0){
    switch(.fmt,
           df  = {return(object@taxa)},
           tbl = {return(tibble::as_tibble(object@taxa, rownames="taxa_id"))})
  }
  
  if(length(object@taxa) == 0 && length(object@comm) != 0){
    switch(.fmt,
           df  = {return(data.frame(comm_id = comm_id(object),
                                    row.names = taxa_id(object)))},
           tbl = {return(comm_id(object, .fmt = "tbl"))})
  }
  
  if(length(object@taxa) != 0 && length(object@comm) != 0){
    
    result <- comm_id(object, .fmt = "tbl") %>%
      left_join(tibble::rownames_to_column(object@taxa, "taxa_id"),  
                by = "taxa_id")
    
    switch(.fmt,
           df  = {return(tibble::column_to_rownames(result, "taxa_id"))},
           tbl = {return(result)})
    
  }
  
  
})

setMethod("taxa", "mgnets", function(object, .fmt = "df", .collapse = FALSE) {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if (.collapse) {
    taxa_collapsed <- object %>% 
      purrr::map(\(x) taxa(x, .fmt = "tbl")) %>% 
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>% 
      purrr::list_rbind()
    return(taxa_collapsed)
  }
  
  if(.fmt == "df") {
    
    if(length(object) == 0) return(data.frame())
    return(sapply(object, taxa, .fmt = "df", simplify = FALSE, USE.NAMES = TRUE))
    
  } else if(.fmt == "tbl") {
    
    if(length(object) == 0) return(tibble::tibble())
    return(sapply(object, taxa, .fmt = "tbl", simplify = FALSE, USE.NAMES = TRUE))
    
  } 
  
})


#' @title Retrieve network graph(s)
#'
#' @description
#' Return the network stored in the \code{netw} slot.
#' For an \code{mgnet}, a single \code{igraph} is returned.
#' For an \code{mgnets}, a named \code{list} of \code{igraph} (one per element).
#' If the \code{taxa} slot is non-empty, its columns are attached as vertex attributes.
#'
#' @param object An \code{mgnet} or \code{mgnets} object.
#' @param selected Logical (default \code{TRUE}). If \code{TRUE}, for each graph
#'   only edges whose \code{link_id} is currently selected (via \code{\link{select_link}})
#'   are kept; otherwise the full graph is returned.
#'
#' @return
#' \itemize{
#'   \item \code{mgnet}: an \code{igraph} object.
#'   \item \code{mgnets}: a named \code{list} of \code{igraph} objects.
#' }
#'
#' @details
#' Requires a non-missing network in the object. When \code{selected = TRUE} and no
#' links are selected, an empty-edge graph is returned for that element.
#'
#' @seealso \code{\link{select_link}}, \code{\link{deselect_link}}, \code{\link{taxa}}
#'
#' @importFrom igraph vertex_attr set_edge_attr is_weighted E membership
#' @aliases netw,mgnet-method netw,mgnets-method
#' @export
setGeneric("netw", function(object, selected = TRUE) standardGeneric("netw"))

setMethod("netw", "mgnet", function(object, selected = TRUE) {
  g <- object@netw
  if (length(g) == 0) return(g)
  
  if (has_metataxa(object)) {
    metataxa <- taxa(object)
    for (vertex_attr in colnames(metataxa)) {
      igraph::vertex_attr(g, vertex_attr) <- metataxa[[vertex_attr]]
    }
  }
  
  if (isTRUE(selected) && are_selected_links(object)) {
    sel_ids <- get_selected_links(object)                  
    all_ids <- igraph::edge_attr(g, "link_id")             
    keep    <- which(all_ids %in% sel_ids)                 
    g <- igraph::subgraph_from_edges(g, eids = keep, delete.vertices = FALSE)
  }
  
  g
})


setMethod("netw", "mgnets", function(object, selected = TRUE) {
  sapply(object@mgnets, function(x) netw(x, selected = selected),
         simplify = FALSE, USE.NAMES = TRUE)
})


# COMMUNITY
#------------------------------------------------------------------------------#
#' Get Community Detection Results
#'
#' Retrieves the community detection results stored in the `comm` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnets`. These results are typically derived from network
#' analysis methods and indicate the grouping or clustering of taxa into communities based on their
#' network interactions.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @return The community detection result object stored in the `comm` slot for an `mgnet` object,
#'         or a list of such objects for an `mgnets` object, each representing the community
#'         detection results of a contained `mgnet` object.
#' @export
#' @name comm
#' @aliases comm,mgnet-method comm,mgnets-method
setGeneric("comm", function(object) standardGeneric("comm"))

setMethod("comm", "mgnet", function(object) {
  object@comm
})

setMethod("comm", "mgnets", function(object) {
  sapply(object@mgnets, function(x) x@comm, simplify = FALSE, USE.NAMES = TRUE)
})


# LINK
#------------------------------------------------------------------------------#
#' Retrieve network links (edges)
#'
#' @description
#' `link()` returns the edge table of the network stored in an `mgnet` or, for an
#' `mgnets`, returns one edge table per contained object (as a named list) or a
#' single stacked tibble with an `mgnet` column.
#'
#' @section Edge table schema:
#' Matches `igraph::as_data_frame(g, what = "edges")`, i.e. at minimum:
#' - `from`, `to` (character vertex ids),
#' - any additional edge attributes (e.g., `weight`, `link_id`, ...).
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param selected Logical scalar. If `TRUE`, use the selected network (if your
#'   class distinguishes selected vs full network). For `mgnet` the default is
#'   `FALSE`. For `mgnets` the default is `TRUE`.
#' @param .fmt For `mgnets` only, output format: one of `"list"` (named list of
#'   tibbles) or `"tbl"` (single tibble with column `mgnet`). Ignored for `mgnet`.
#'
#' @return
#' - For an `mgnet`: a tibble/data frame of edges.
#' - For an `mgnets` with `.fmt = "list"`: a **named list** of edge tibbles.
#' - For an `mgnets` with `.fmt = "tbl"`: a single tibble with an `mgnet` column.
#'
#' @name link
#' @export
setGeneric("link", function(object, selected = FALSE, .fmt = c("list", "tbl"))
  standardGeneric("link"))

#' @rdname link
#' @export
setMethod("link", "mgnet", function(object, selected = FALSE, .fmt = c("list", "tbl")) {
  # validate `selected`
  if (!(is.logical(selected) && length(selected) == 1L && !is.na(selected))) {
    cli::cli_abort(c(
      "x" = "{.arg selected} must be a single {.cls logical} (TRUE/FALSE) and not NA.",
      "i" = "Got class {.cls {class(selected)[1]}} with length {length(selected)}."
    ))
  }
  # network required
  if (miss_netw(object)) cli::cli_abort("No network available in this {.cls mgnet} object.")
  
  g <- netw(object, selected = selected)
  edges <- igraph::as_data_frame(g, what = "edges")
  tibble::as_tibble(edges)
})

#' @rdname link
#' @export
setMethod("link", "mgnets", function(object, selected = FALSE, .fmt = c("list", "tbl")) {
  # validate `selected`
  if (!(is.logical(selected) && length(selected) == 1L && !is.na(selected))) {
    cli::cli_abort(c(
      "x" = "{.arg selected} must be a single {.cls logical} (TRUE/FALSE) and not NA.",
      "i" = "Got class {.cls {class(selected)[1]}} with length {length(selected)}."
    ))
  }
  .fmt <- match.arg(.fmt, c("list", "tbl"))
  
  if (miss_netw(object, "any")) {
    cli::cli_abort("No network available in at least one element of this {.cls mgnets}.")
  }
  
  res <- sapply(object@mgnets, function(x) link(x, selected = selected),
                simplify = FALSE, USE.NAMES = TRUE)
  if (.fmt == "tbl") res <- purrr::list_rbind(res, names_to = "mgnet")
  
  res
})




