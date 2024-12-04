# ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Abundance Data Based on Taxonomic Information
#'
#' @description
#' Retrieves and aggregates abundance data for an `mgnet` or `mgnetList` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnetList` object containing abundance and informational data.
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original abundance data without any aggregation.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .fun A function to specify how the abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnetList` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name abun
#' @aliases abun,mgnet-method abun,mgnetList-method
setGeneric("abun", function(object, .var = NULL, .fmt = "mat", .fun = sum) standardGeneric("abun"))

setMethod("abun", "mgnet", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  
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
    
    taxonomic_info <- taxa(object, .fmt = "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
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

setMethod("abun", "mgnetList", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@mgnets, function(x) abun(object = x, .var = .var, 
                                                        .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# RELATIVE ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Relative Abundance Data Based on Taxonomic Information
#'
#' @description
#' Retrieves and aggregates relative abundance data for an `mgnet` or `mgnetList` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original relative abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain relative abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnetList` object containing relative abundance and informational data.
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original relative abundance data without any aggregation.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the relative abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns relative abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnetList` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name rela
#' @aliases rela,mgnet-method rela,mgnetList-method
setGeneric("rela", function(object, .var = NULL, .fmt = "mat", .fun = sum) standardGeneric("rela"))

setMethod("rela", "mgnet", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  
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
    
    taxonomic_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
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

setMethod("rela", "mgnetList", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- sapply(object@mgnets, function(x) rela(object = x, .var = .var, 
                                                            .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# NORMALIZED ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Normalized Abundance Data Based on Taxonomic Information
#'
#' @description
#' Retrieves and aggregates normalized abundance data for an `mgnet` or `mgnetList` object based on
#' a specified column in the `taxa` data frame. This function allows for custom aggregation 
#' methods such as sum or mean. If no column is specified, it returns the original normalized abundance data 
#' in the selected format. Although `taxa` can contain various types of information, it is 
#' commonly used to obtain normalized abundances at different taxonomic ranks or other categorical variables.
#'
#' @param object An `mgnet` or `mgnetList` object containing normalized abundance and informational data.
#' @param .var A character string specifying the column name from the `taxa` slot to be used 
#'        for grouping and aggregation. If this parameter is omitted, the function will return 
#'        the original normalized abundance data without any aggregation.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are:
#'        - "mat": returns a matrix
#'        - "df": returns a data.frame
#'        - "tbl": returns a tibble, where for `mgnet` objects, the row names of the normalized abundance 
#'          matrix are moved into a new column named `sample_id`.
#'        The default format is "mat".
#' @param .fun A function to specify how the relative abundance data should be aggregated. 
#'        This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
#'        The default aggregation function is sum.
#'
#' @return Depending on the '.fmt' parameter:
#'         - For an `mgnet` object, returns normalized abundance data formatted as specified, aggregated
#'           based on the provided column if specified.
#'         - For an `mgnetList` object, returns a list of such formatted data from each `mgnet`
#'           object within the list, enabling batch processing and analysis of multiple datasets.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise arrange
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @name norm
#' @aliases norm,mgnet-method norm,mgnetList-method
setGeneric("norm", function(object, .var = NULL, .fmt = "mat", .fun = sum) standardGeneric("norm"))

setMethod("norm", "mgnet", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  
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
    
    taxonomic_info <- taxa(object, "tbl") %>%
      dplyr::select(taxa_id, !!rlang::sym(.var))
    
    # Aggregate data based on .var
    aggregated_data <- data_frame %>%
      dplyr::left_join(taxonomic_info, by = "taxa_id") %>%
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

setMethod("norm", "mgnetList", function(object, .var = NULL, .fmt = "mat", .fun = sum) {
  
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
#' or for each `mgnet` object within an `mgnetList`, with the option to format the output as
#' a `data.frame`, `tibble`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, the row names of
#'          the abundance matrix are converted into a new column named `sample_id`, ensuring alignment
#'          with the reserved keyword in `mgnet-class`.
#'          
#' @return The content of the `meta` slot for `mgnet` or a list of such contents for `mgnetList`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name meta
#' @aliases meta,mgnet-method meta,mgnetList-method
setGeneric("meta", function(object, .fmt = "df") standardGeneric("meta"))

setMethod("meta", "mgnet", function(object, .fmt = "df") {
  
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

setMethod("meta", "mgnetList", function(object, .fmt = "df") {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
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
#' or for each `mgnet` object within an `mgnetList`, with the option to format the output as
#' a `data.frame`, `tibble`, or a combined `tibble` for multiple mgnet objects.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, the the row names of
#'          the abundance matrix are converted into a new column named `taxa_id`, ensuring alignment
#'          with the reserved keyword in `mgnet-class`.
#'          
#' @return The content of the `taxa` slot for `mgnet` or a list of such contents for `mgnetList`.
#' 
#' @importFrom purrr map imap list_rbind
#' @importFrom tibble tibble
#' @export
#' @name taxa
#' @aliases taxa,mgnet-method taxa,mgnetList-method
setGeneric("taxa", function(object, .fmt = "df") standardGeneric("taxa"))

setMethod("taxa", "mgnet", function(object, .fmt = "df") {
  
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
           df  = {return(comm_id(object, .fmt = "df"))},
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

setMethod("taxa", "mgnetList", function(object, .fmt = "df") {
  
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if(.fmt == "df") {
    
    if(length(object) == 0) return(data.frame())
    return(sapply(object, taxa, .fmt = "df", simplify = FALSE, USE.NAMES = TRUE))
    
  } else if(.fmt == "tbl") {
    
    if(length(object) == 0) return(tibble::tibble())
    return(sapply(object, taxa, .fmt = "tbl", simplify = FALSE, USE.NAMES = TRUE))
    
  } 
  
})


# NETWORK
#------------------------------------------------------------------------------#
#' Retrieve Network Graph from mgnet Objects
#'
#' Retrieves the network graph stored in the `network` slot of an `mgnet` object
#' or from each `mgnet` object within an `mgnetList`. This network graph is typically an `igraph` object
#' representing the interactions among different taxa.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param add_vertex_attr A logical indicating whether to attach all available taxa information as attributes of the vertices. Defaults to FALSE.
#' @return An `igraph` object containing the network from an `mgnet` object, or a list of `igraph` objects
#'         from an `mgnetList`, each representing the network graph of a contained `mgnet` object.
#' @export
#' @importFrom igraph vertex_attr set_edge_attr is_weighted E membership
#' @aliases netw,mgnet-method netw,mgnetList-method
setGeneric("netw", function(object, add_vertex_attr = FALSE) standardGeneric("netw"))

setMethod("netw", "mgnet", function(object, add_vertex_attr = FALSE) {

  g <- object@netw
  if (length(g) == 0) return(g)

  if (add_vertex_attr) {

    # Add vertex attributes
    if (length(object@taxa) != 0) {
      for (vertex_attr in names(object@taxa)) {
        igraph::vertex_attr(g, vertex_attr) <- object@taxa[[vertex_attr]]
      }
    }
    
    if (length(object@comm) != 0){
      igraph::vertex_attr(g, "comm_id") <- igraph::membership(object@comm)
    }
  }

  return(g)
})


setMethod("netw", "mgnetList", function(object, add_vertex_attr) {
  sapply(object@mgnets, function(x) netw(x, add_vertex_attr), simplify = FALSE, USE.NAMES = TRUE)
})


# LINK
#------------------------------------------------------------------------------#
#' Retrieve Edge List with Metadata from mgnet Object(s)
#'
#' Extracts an edge list from the network slot of an `mgnet` or `mgnetList` object, including edge attributes 
#' and associated taxa metadata for both vertices in each link.
#'
#' @param object An `mgnet` or `mgnetList` object.
#'               For `mgnet`, the method will extract the network's edge list and merge taxa metadata 
#'               for both source and target nodes. For `mgnetList` the method will be applied on each element.
#' @param .suffix A character vector of length 2 providing suffixes to append to the taxa 
#'        metadata columns to distinguish the nodes connected from an edge. 
#'        Default is c("_1", "_2"). The values must be distinct to prevent column name overlap.
#'
#' @details
#' This method leverages the network structure within `mgnet` objects to generate a detailed edge tibble that includes:
#' - `Node Identifiers`: Each identifier for the nodes connected by an edge begins with the string "taxa" and is 
#'   appended with `.suffix` to indicate the source and target nodes, respectively (e.g., `taxa_id_1`, `taxa_id_2`).
#' - `Link Attributes`: Includes all available attributes associated with the links, such as weight.
#' - `Metadata`: Metadata from both source and target taxa are included, with column names appended 
#'   with `_1` and `_2` as suffixes. These suffixes help distinguish between the source and target taxa metadata.
#'
#' For `mgnetList` objects, this transformation is applied individually to each `mgnet` object.
#'
#' @importFrom igraph E is_weighted as_data_frame
#' @importFrom dplyr left_join
#' @importFrom purrr map imap list_rbind
#' @export
#' @name link
#' @aliases link,mgnet-method link,mgnetList-method
setGeneric("link", function(object, .suffix = c("_1", "_2")) standardGeneric("link"))

setMethod("link", "mgnet", function(object, .suffix = c("_1", "_2")) {
  
  # Ensure the network is available
  if (miss_slot(object, "netw")) stop("Error: No network available.")
  
  # Check .suffix
  if(length(.suffix) != 2 || !is.character(.suffix) || .suffix[1] == .suffix[2]){
    stop(".suffix must be a character vector of length 2 with distinct values")
  }
  
  # Check if is requirred the edge filter
  selected_links <- get_selected_links(object)
  if (!is.null(selected_links)) {
    netw0 <- netw(object)
    netw(object) <- igraph::subgraph_from_edges(graph = netw(object),
                                                eids = get_selected_links(object),
                                                delete.vertices = FALSE)
  }
  
  # Extract edge list with weights
  net <- netw(object)
  edges_df <- igraph::as_data_frame(net, what = "edges")
  edges_df <- as_tibble(edges_df)
  colnames(edges_df)[1:2] <- paste0("taxa_id", .suffix)
  
  new_taxa_info_cols <- paste0(taxa_vars(object), .suffix[1])
  new_taxa_info_cols <- c(new_taxa_info_cols, paste0(taxa_vars(object), .suffix[2]))
  
  # Check for potential column name overlaps
  if (any(new_taxa_info_cols %in% names(edges_df)[-c(1,2)])) {
    overlapping_columns <- new_taxa_info_cols[new_taxa_info_cols %in% names(edges_df)]
    stop("Overlap detected in column names between nodes and edges: ", paste(overlapping_columns, collapse = ", "))
  }

  # First, ensure you get the network correctly and extract attribute names
  edges__attr_names <- names(igraph::edge_attr(netw(object)))
  
  if (!is.null(edges__attr_names) && length(edges__attr_names) > 0) {
    
    # Create a regex pattern to match any of the suffixes at the end of strings
    pattern <- stringr::str_c(.suffix, collapse = "|", suffix = "$")
    
    # Use str_detect to find elements ending with any suffix
    invalid_names <- edges__attr_names[stringr::str_detect(edges__attr_names, pattern)]
    if (length(invalid_names) > 0) {
      # Construct a concise error message that lists the problematic attribute names and the disallowed suffixes
      suffix_list <- paste(.suffix, collapse = ", ")
      error_message <- sprintf(
        "Edge attributes '%s' end with disallowed suffixes (%s). " +
          "Suffixes must differ between nodes and edges to discern them. Please adjust the .suffix values or rename the edge or vertex attributes accordingly.",
        paste(invalid_names, collapse = ", "), suffix_list
      )
      stop(error_message)
    }
  }
  
  # Merge with taxa metadata for source nodes
  tbl_1 <- gather_taxa(object)
  colnames(tbl_1) <- paste0(colnames(tbl_1), .suffix[1])
  edges_df <- dplyr::left_join(edges_df, tbl_1, by = paste0("taxa_id", .suffix[1]))
  
  # Merge with taxa metadata for target nodes
  tbl_2 <- gather_taxa(object)
  colnames(tbl_2) <- paste0(colnames(tbl_2), .suffix[2])
  edges_df <- dplyr::left_join(edges_df, tbl_2, by = paste0("taxa_id", .suffix[2]))
  
  return(edges_df)
  
})

setMethod("link", "mgnetList", function(object, .suffix = c("_1", "_2")){
  sapply(object, \(x) link(x, .suffix), simplify = FALSE, USE.NAMES = TRUE)
})


# COMMUNITY
#------------------------------------------------------------------------------#
#' Get Community Detection Results
#'
#' Retrieves the community detection results stored in the `comm` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`. These results are typically derived from network
#' analysis methods and indicate the grouping or clustering of taxa into communities based on their
#' network interactions.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return The community detection result object stored in the `comm` slot for an `mgnet` object,
#'         or a list of such objects for an `mgnetList` object, each representing the community
#'         detection results of a contained `mgnet` object.
#' @export
#' @name comm
#' @aliases comm,mgnet-method comm,mgnetList-method
setGeneric("comm", function(object) standardGeneric("comm"))

setMethod("comm", "mgnet", function(object) {
  object@comm
})

setMethod("comm", "mgnetList", function(object) {
  sapply(object@mgnets, function(x) x@comm, simplify = FALSE, USE.NAMES = TRUE)
})

# META_VARS
#------------------------------------------------------------------------------#
#' Get Sample Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `meta` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of sample metadata variables 
#'        in each `mgnet` and \code{"unique"} for a single array with the unique elemnts. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#'         
#' @examples
#' data(mg, package = "mgnet")
#' meta_vars(mg)  
#'
#' data(mgl, package = "mgnet")
#' meta_vars(mgl, .fmt = "list)  
#' meta_vars(mgl, .fmt = "unique) 
#' 
#' @export
#' @name meta_vars
#' @aliases meta_vars,mgnet-method meta_vars,mgnetList-method
setGeneric("meta_vars", function(object, .fmt = "list") standardGeneric("meta_vars"))

setMethod("meta_vars", "mgnet", function(object, .fmt = "list") {
  
  if(miss_sample(object)){
    return(character(length=0))
  } else if(length(object@meta) == 0){
    return("sample_id")
  } else {
    return(colnames(meta(object, "tbl")))
  }
  
})

setMethod("meta_vars", "mgnetList", function(object, .fmt = "list") {
  
  .fmt <- match.arg(.fmt, c("list", "unique"))
  
  vars <- sapply(object, meta_vars, simplify = FALSE, USE.NAMES = TRUE)
  vars <- sapply(vars, \(x){if(length(x)!=0) c("mgnet", x) else x}, 
                 simplify = FALSE, USE.NAMES = TRUE)
  
  if(.fmt == "list"){
    return(vars)
  } else {
    return(unique(unlist(vars)))
  }
  
})

# TAXA_VARS
#------------------------------------------------------------------------------#
#' Get Taxa Metadata Variables
#'
#' Retrieves the names of metadata variables available in the `taxa` slot of an `mgnet` object,
#' or for each `mgnet` object within an `mgnetList`.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt Character; specifies the output format when checking an `mgnetList`. 
#'        Accepted values are \code{"list"} for a list of sample metadata variables 
#'        in each `mgnet` and \code{"unique"} for a single array with the unique elemnts. 
#'        Default is \code{"list"}.
#'        
#' @return For an `mgnet` object, a character vector of metadata variable names.
#'         For an `mgnetList` object, a named list of character vectors, with each list item representing 
#'         the metadata variable names in the corresponding `mgnet` objects.
#'         
#' @examples
#' data(mg, package = "mgnet")
#' taxa_vars(mg)  
#'
#' data(mgl, package = "mgnet")
#' taxa_vars(mgl, .fmt = "list)  
#' taxa_vars(mgl, .fmt = "unique) 
#' 
#' @export
#' @name taxa_vars
#' @aliases taxa_vars,mgnet-method taxa_vars,mgnetList-method
setGeneric("taxa_vars", function(object, .fmt = "list") standardGeneric("taxa_vars"))

setMethod("taxa_vars", "mgnet", function(object, .fmt = "list") {
  
  if(!has_sample(object)){
    return(character(length=0))
  } else if(miss_metataxa(object)){
    return("taxa_id")
  } else {
    return(colnames(taxa(object, "tbl")))
  }
  
})

setMethod("taxa_vars", "mgnetList", function(object, .fmt = "list") {
  
  .fmt <- match.arg(.fmt, c("list", "unique"))
  
  vars <- sapply(object, taxa_vars, simplify = FALSE, USE.NAMES = TRUE)
  vars <- sapply(vars, \(x){if(length(x)!=0) c("mgnet", x) else x}, 
                 simplify = FALSE, USE.NAMES = TRUE)
  
  if(.fmt == "list"){
    return(vars)
  } else {
    return(unique(unlist(vars)))
  }
  
})

