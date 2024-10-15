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
  result <- sapply(object@mgnets, function(x) norm(object = x, .var = .var, 
                                                             .fmt = .fmt, .fun = .fun),
                   simplify = FALSE, USE.NAMES = TRUE)
  return(result)
})


# INFO_SAMPLE
#------------------------------------------------------------------------------#
#' Get Sample Information
#'
#' Retrieves the sample information stored in the `meta` slot of an `mgnet` object
#' or for each `mgnet` object within an `mgnetList`, with the option to format the output as
#' a `data.frame`, `tibble`, or a combined `tibble` for multiple mgnet objects.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, if present, the row names of
#'          the abundance matrix are converted into a new column named `sample_id`, ensuring alignment
#'          with the reserved keyword in `mgnet`.
#'        - "list_df": When working with `mgnetList`, returns a list of `data.frame` objects, 
#'          each corresponding to the `meta` of an individual `mgnet` object.
#'        - "list_tbl": Similar to "list_df", but each entry in the list is a `tibble`.
#'        - "tbl" (for `mgnetList` only): Concatenates the results of all `mgnet` objects into a single 
#'          `tibble`, with an additional column identifying the source `mgnet` object. The column is named 
#'          `mgnet`, which is a reserved keyword of the `mgnet` class and contains the name of each `mgnet` object.
#'          
#' @return The content of the `meta` slot for `mgnet` or a list of such contents for `mgnetList`.
#' 
#' @importFrom purrr map imap list_rbind
#' @export
#' @name meta
#' @aliases meta,mgnet-method meta,mgnetList-method
setGeneric("meta", function(object, .fmt) standardGeneric("meta"))

setMethod("meta", "mgnet", function(object, .fmt) {
  
  if(missing(.fmt)) .fmt<- "df"
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

setMethod("meta", "mgnetList", function(object, .fmt) {
  if(missing(.fmt)) .fmt <- "list_df"
  .fmt <- match.arg(.fmt, c("list_df", "list_tbl", "tbl"))
  
  if(.fmt == "list_df") {
    
    return(purrr::map(object, function(x) meta(x, .fmt = "df")))
    
  } else if(.fmt == "list_tbl") {
    
    return(purrr::map(object, function(x) meta(x, .fmt = "tbl")))
    
  } else {
    
    purrr::map(object, function(x) meta(x, .fmt = "tbl")) %>%
      purrr::imap(\(x,y) mutate(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      return()
    
  }
  
})


# INFO_TAXA
#------------------------------------------------------------------------------#
#' Get Taxa Information
#'
#' Retrieves the taxa information stored in the `taxa` slot of an `mgnet` object
#' or for each `mgnet` object within an `mgnetList`, with the option to format the output as
#' a `data.frame`, `tibble`, or a combined `tibble` for multiple mgnet objects.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result.
#'        Possible choices are:
#'        - "df": Returns the output as a `data.frame`.
#'        - "tbl": Returns the output as a `tibble`. For `mgnet` objects, if present, the row names of
#'          the abundance matrix are converted into a new column named `taxa_id`, ensuring alignment
#'          with the reserved keyword in `mgnet`.
#'        - "list_df": When working with `mgnetList`, returns a list of `data.frame` objects, 
#'          each corresponding to the `taxa` of an individual `mgnet` object.
#'        - "list_tbl": Similar to "list_df", but each entry in the list is a `tibble`.
#'        - "tbl" (for `mgnetList` only): Concatenates the results of all `mgnet` objects into a single 
#'          `tibble`, with an additional column identifying the source `mgnet` object. The column is named 
#'          `mgnet`, which is a reserved keyword of the `mgnet` class and contains the name of each `mgnet` object.
#'          
#' @return The content of the `taxa` slot for `mgnet` or a list of such contents for `mgnetList`.
#' 
#' @importFrom purrr map imap list_rbind
#' @export
#' @name taxa
#' @aliases taxa,mgnet-method taxa,mgnetList-method
setGeneric("taxa", function(object, .fmt) standardGeneric("taxa"))

setMethod("taxa", "mgnet", function(object, .fmt) {
  
  if(missing(.fmt)) .fmt <- "df"
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

setMethod("taxa", "mgnetList", function(object, .fmt) {
  
  if(missing(.fmt)) .fmt <- "list_df"
  .fmt <- match.arg(.fmt, c("list_df", "list_tbl", "tbl"))
  
  if(.fmt == "list_df") {
    
    return(purrr::map(object, function(x) taxa(x, .fmt = "df")))
    
  } else if(.fmt == "list_tbl") {
    
    return(purrr::map(object, function(x) taxa(x, .fmt = "tbl")))
    
  } else {
    
    purrr::map(object, function(x) taxa(x, .fmt = "tbl")) %>%
      purrr::imap(\(x,y) mutate(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind() %>%
      return()
    
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
#' @importFrom igraph vertex_attr set_edge_attr is_weighted E
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
  }

  return(g)
})


setMethod("netw", "mgnetList", function(object, add_vertex_attr) {
  sapply(object@mgnets, function(x) netw(x, add_vertex_attr), simplify = FALSE, USE.NAMES = TRUE)
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
