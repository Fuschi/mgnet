# ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Abundance Data at Specified Taxonomic Rank
#'
#' @description
#' Retrieves abundance data for an `mgnet` or `mgnetList` object, optionally 
#' aggregated at a specified taxonomic rank, and allows the output to be formatted 
#' as a matrix, data frame, or tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param rank A character string specifying the taxonomic rank of interest for data aggregation.
#'        If not provided or if the rank is "missing", the function returns the original abundance 
#'        data without aggregation. 
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "mat" for matrix, "df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`.
#'
#' @return For an `mgnet` object, abundance data formatted as specified by the `.fmt` parameter,
#'         aggregated at the specified taxonomic rank if provided. For an `mgnetList` object, 
#'         a list of such formatted data from each `mgnet` object within the list, allowing 
#'         for batch processing and analysis of multiple datasets simultaneously.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column tibble
#' @importFrom tidyr pivot_longer pivot_wider
#' @name abundance
#' @aliases abundance,mgnet-method abundance,mgnetList-method
setGeneric("abundance", function(object, rank = "missing", .fmt = "mat") standardGeneric("abundance"))

setMethod("abundance", "mgnet", function(object, rank, .fmt) {
  
  # Checks
  if (missing(rank)) rank <- "missing"
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  if( length(object@abundance) == 0 ){
    switch(.fmt,
           mat = {return(matrix(nrow=0, ncol=0))},
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  } else if(rank=="missing") {
    switch(.fmt,
           mat = {return(object@abundance)},
           df  = {return(as.data.frame(object@abundance))},
           tbl = {return(tibble::as_tibble(object@abundance, rownames="sample_id"))})
  }  else {
    if(!rank %in% colnames(object@lineage)) {
      stop("rank must be present in the object. Available ranks are: ",
           paste(toString(colnames(object@lineage)), collapse=", "))}
    
    data_frame <- object@abundance %>% tibble::as_tibble(rownames="sample_id") %>%
      tidyr::pivot_longer(cols=-sample_id, names_to="taxa_id", values_to="abundance")
      
    lineage_info <- object@lineage %>%
      tibble::as_tibble(rownames="taxa_id") %>%
      dplyr::select(taxa_id, !!rlang::sym(rank))
      
    aggregated_data <- data_frame %>%
      dplyr::left_join(lineage_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(rank)) %>%
      dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(rank), values_from = abundance) %>%
      dplyr::arrange(sample_id = sample_id(object))
      
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

setMethod("abundance", "mgnetList", function(object, rank = "missing", .fmt = "mat") {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- lapply(object@mgnets, function(x) abundance(x, rank, .fmt))
  return(result)
})


# REL_ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Relative Abundance Data at Specified Taxonomic Rank
#'
#' @description
#' Retrieves relative abundance data for an `mgnet` or `mgnetList` object, optionally 
#' aggregated at a specified taxonomic rank, and allows the output to be formatted 
#' as a matrix, data frame, or tibble.
#' 
#' @details
#' The relative abundance matrix is obtained by dividing each sample's abundance values by 
#' the corresponding sample sum present in the \code{sample_sum} column of the \code{info_sample}
#' data.frame associated with the mgnet object. The abundance matrices at the chosen rank are
#' retrieved using the `\link{abundance}` method.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param rank A character string specifying the taxonomic rank of interest for data aggregation.
#'        If not provided or if the rank is "missing", the function returns the original abundance 
#'        data without aggregation. 
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "mat" for matrix, "df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`.
#'
#' @return For an `mgnet` object, abundance data formatted as specified by the `.fmt` parameter,
#'         aggregated at the specified taxonomic rank if provided. For an `mgnetList` object, 
#'         a list of such formatted data from each `mgnet` object within the list, allowing 
#'         for batch processing and analysis of multiple datasets simultaneously.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider
#' @name rel_abundance
#' @aliases rel_abundance,mgnet-method rel_abundance,mgnetList-method
setGeneric("rel_abundance", function(object, rank = "missing", .fmt = "mat") standardGeneric("rel_abundance"))

setMethod("rel_abundance", "mgnet", function(object, rank, .fmt) {
  
  # Checks
  if (missing(rank)) rank <- "missing"
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  if( length(object@rel_abundance) == 0 ){
    switch(.fmt,
           mat = {return(matrix(nrow=0, ncol=0))},
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  } else if(rank=="missing"){
    switch(.fmt,
           mat = {return(object@rel_abundance)},
           df  = {return(as.data.frame(object@rel_abundance))},
           tbl = {return(tibble::as_tibble(object@rel_abundance, rownames="sample_id"))})
  } else {
    if(!rank %in% colnames(object@lineage)) {
      stop("rank must be present in the object. Available ranks are: ",
           paste(toString(colnames(object@lineage)), collapse=", "))}
    
    data_frame <- object@rel_abundance %>% tibble::as_tibble(rownames="sample_id") %>%
      tidyr::pivot_longer(cols=-sample_id, names_to="taxa_id", values_to="rel_abundance")
    
    lineage_info <- object@lineage %>%
      tibble::as_tibble(rownames="taxa_id") %>%
      dplyr::select(taxa_id, !!rlang::sym(rank))
    
    aggregated_data <- data_frame %>%
      dplyr::left_join(lineage_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(rank)) %>%
      dplyr::summarise(rel_abundance = sum(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(rank), values_from = rel_abundance) %>%
      dplyr::arrange(sample_id = sample_id(object))
    
    result <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
    
    switch(.fmt,
           mat = {return(result)},
           df  = {return(as.data.frame(result))},
           tbl = {return(tibble::as_tibble(result, rownames="sample_id"))}
    )
  }
})

setMethod("rel_abundance", "mgnetList", function(object, rank = "missing", .fmt = "mat") {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- lapply(object@mgnets, function(x) rel_abundance(x, rank, .fmt))
  return(result)
})


# NORM_ABUNDANCE
#------------------------------------------------------------------------------#
#' Retrieve Log_Abundance Data at Specified Taxonomic Rank
#'
#' @description
#' Retrieves log-ratio transformed abundance data for an `mgnet` or `mgnetList` object, 
#' optionally aggregated at a specified taxonomic rank, and allows the output to 
#' be formatted as a matrix, data frame, or tibble.
#' 
#' @details
#' This function returns the log-ratio transformed abundance data, which is useful 
#' for compositional data analysis.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param rank A character string specifying the taxonomic rank of interest for data aggregation.
#'        If not provided or if the rank is "missing", the function returns the original abundance 
#'        data without aggregation. 
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "mat" for matrix, "df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`.
#'
#' @return For an `mgnet` object, abundance data formatted as specified by the `.fmt` parameter,
#'         aggregated at the specified taxonomic rank if provided. For an `mgnetList` object, 
#'         a list of such formatted data from each `mgnet` object within the list, allowing 
#'         for batch processing and analysis of multiple datasets simultaneously.
#'
#' @export
#' @importFrom dplyr %>% select left_join group_by summarise
#' @importFrom rlang sym
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider
#' @name norm_abundance
#' @aliases norm_abundance,mgnet-method norm_abundance,mgnetList-method
setGeneric("norm_abundance", function(object, rank = "missing", .fmt = "mat") standardGeneric("norm_abundance"))

setMethod("norm_abundance", "mgnet", function(object, rank, .fmt) {
  
  # Checks
  if (missing(rank))rank <- "missing"
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  
  if( length(object@norm_abundance) == 0 ){
    
    switch(.fmt,
           mat = {return(matrix(nrow=0, ncol=0))},
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
    
  } else if(rank=="missing"){
    
    result <- object@norm_abundance
    switch(.fmt,
           mat = {return(result)},
           df  = {return(as.data.frame(result))},
           tbl = {return(tibble::as_tibble(result, rownames="sample_id"))})
    
  } else {
    
    if( !rank %in% colnames(object@lineage) ) {
      stop("specified rank must be present in the object. Available ranks are: ",
           paste(toString(colnames(object@lineage)), collapse=", "))}
    
    data_frame <- object@norm_abundance %>% tibble::as_tibble(rownames="sample_id") %>%
      tidyr::pivot_longer(cols=-sample_id, names_to="taxa_id", values_to="norm_abundance")
    
    lineage_info <- object@lineage %>%
      tibble::as_tibble(rownames="taxa_id") %>%
      dplyr::select(taxa_id, !!rlang::sym(rank))
    
    aggregated_data <- data_frame %>%
      dplyr::left_join(lineage_info, by = "taxa_id") %>%
      dplyr::group_by(sample_id, !!rlang::sym(rank)) %>%
      dplyr::summarise(norm_abundance = sum(norm_abundance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(rank), values_from = norm_abundance) %>%
      dplyr::arrange(sample_id = sample_id(object))
    
    result <- aggregated_data %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()

    switch(.fmt,
           mat = {return(result)},
           df  = {return(as.data.frame(result))},
           tbl = {return(tibble::as_tibble(result, rownames="sample_id"))}
    )
  }
})

setMethod("norm_abundance", "mgnetList", function(object, rank = "missing", .fmt = "mat") {
  .fmt <- match.arg(.fmt, c("mat", "df", "tbl"))
  result <- lapply(object@mgnets, function(x) norm_abundance(x, rank, .fmt))
  return(result)
})



# INFO_SAMPLE
#------------------------------------------------------------------------------#
#' Get Sample Information
#'
#' Retrieves the sample information stored in the `info_sample` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`, with the option to format the output as
#' a data.frame or tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`.
#' @return The content of the `info_sample` slot for `mgnet` or a list of such contents for `mgnetList`.
#' @export
#' @name info_sample
#' @aliases info_sample,mgnet-method info_sample,mgnetList-method
setGeneric("info_sample", function(object, .fmt = "df") standardGeneric("info_sample"))

setMethod("info_sample", "mgnet", function(object, .fmt = "df") {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if( length(object@info_sample) == 0 ){
    
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  } else {
    
    switch(.fmt,
           df  = {return(object@info_sample)},
           tbl = {return(tibble::as_tibble(object@info_sample, rownames="sample_id"))})
    
  }
})

setMethod("info_sample", "mgnetList", function(object, .fmt = "df") {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  sapply(object@mgnets, function(x) info_sample(x, .fmt),
         simplify = FALSE, USE.NAMES = TRUE)
})


# LINEAGE
#------------------------------------------------------------------------------#
#' Get Lineage Information
#'
#' Retrieves the taxonomic lineage information stored in the `lineage` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`, with the option to format the output as
#' matrix, data.frame or tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "mat" for matrix, df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`. Default is "mat".
#' @return The content of the `lineage` slot for `mgnet` or a list of such contents for `mgnetList`.
#' @export
#' @name lineage
#' @aliases lineage,mgnet-method lineage,mgnetList-method
setGeneric("lineage", function(object, .fmt = "mat") standardGeneric("lineage"))

setMethod("lineage", "mgnet", function(object, .fmt = "mat") {
  
  .fmt <- match.arg(.fmt, c("mat","df", "tbl"))
  
  if( length(object@lineage) == 0 ){
    
    switch(.fmt,
           mat = {return(matrix(character(0),nrow=0,ncol=0))},
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  } else {
    
    switch(.fmt,
           mat  = {return(object@lineage)},
           df  = {return(data.frame(object@lineage))},
           tbl = {return(tibble::as_tibble(object@lineage, rownames="taxa_id"))})
    
  }
  
})

setMethod("lineage", "mgnetList", function(object, .fmt = "mat") {
  .fmt <- match.arg(.fmt, c("mat","df", "tbl"))
  sapply(object@mgnets, function(x, .fmt) lineage(x,.fmt),
         simplify = FALSE, USE.NAMES = TRUE)
})



# INFO_TAXA
#------------------------------------------------------------------------------#
#' Get Taxa Information
#' 
#' Retrieves the taxa information stored in the `info_taxa` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`, with the option to format the output as
#' data.frame or tibble.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param .fmt A character string specifying the output format of the result. 
#'        Possible choices are "df" for data.frame, and "tbl" for tibble. 
#'        When ".fmt" is set to "tbl", for `mgnet` objects, the row names of the abundance matrix 
#'        are moved into a new column named `sample_id`, aligning with the reserved keyword in 
#'        `mgnet`.
#' @return The content of the `info_taxa` slot for `mgnet` or a list of such contents for `mgnetList`.
#' @export
#' @name info_taxa
#' @aliases info_taxa,mgnet-method info_taxa,mgnetList-method
setGeneric("info_taxa", function(object, .fmt = "df") standardGeneric("info_taxa"))

setMethod("info_taxa", "mgnet", function(object, .fmt = "df") {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  
  if( length(object@info_taxa) == 0 ){
    
    switch(.fmt,
           df  = {return(data.frame())},
           tbl = {return(tibble::tibble())})
  } else {
    
    switch(.fmt,
           df  = {return(object@info_taxa)},
           tbl = {return(tibble::as_tibble(object@info_taxa, rownames="taxa_id"))})
    
  }
})

setMethod("info_taxa", "mgnetList", function(object, .fmt = "df") {
  .fmt <- match.arg(.fmt, c("df", "tbl"))
  sapply(object@mgnets, function(x) info_taxa(x, .fmt),
         simplify = FALSE, USE.NAMES = TRUE)
})


# NETWORK
#------------------------------------------------------------------------------#
#' Get Network Graph
#'
#' Retrieves the network graph information stored in the `network` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`. This network graph is typically an `igraph` object
#' representing the relationships or interactions among different taxa.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return The `igraph` object stored in the `network` slot for an `mgnet` object, or a list of `igraph` objects
#'         for an `mgnetList` object, each representing the network graph of a contained `mgnet` object.
#' @export
#' @name network
#' @aliases network,mgnet-method network,mgnetList-method
setGeneric("network", function(object) standardGeneric("network"))

setMethod("network", "mgnet", function(object) {
  object@network
})

setMethod("network", "mgnetList", function(object) {
  sapply(object@mgnets, function(x) x@network, simplify = FALSE, USE.NAMES = TRUE)
})


# COMMUNITY
#------------------------------------------------------------------------------#
#' Get Community Detection Results
#'
#' Retrieves the community detection results stored in the `community` slot of an `mgnet` object
#' or each `mgnet` object within an `mgnetList`. These results are typically derived from network
#' analysis methods and indicate the grouping or clustering of taxa into communities based on their
#' network interactions.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @return The community detection result object stored in the `community` slot for an `mgnet` object,
#'         or a list of such objects for an `mgnetList` object, each representing the community
#'         detection results of a contained `mgnet` object.
#' @export
#' @name community
#' @aliases community,mgnet-method community,mgnetList-method
setGeneric("community", function(object) standardGeneric("community"))

setMethod("community", "mgnet", function(object) {
  object@community
})

setMethod("community", "mgnetList", function(object) {
  sapply(object@mgnets, function(x) x@community, simplify = FALSE, USE.NAMES = TRUE)
})


# GETTERS DOCUMENTATION (IT DOES NOT WORKS)
#------------------------------------------------------------------------------#
#' Getter Functions for `mgnet` and `mgnetList` Objects
#'
#' @description
#' The `mgnet` package provides a suite of getter functions designed to extract various 
#' types of data and information from `mgnet` and `mgnetList` objects. These functions 
#' allow users to efficiently access specific components of their metagenomic network 
#' analysis results, facilitating further analysis, visualization, or summary.
#'
#' @details
#' The getter functions cover a range of data types within the `mgnet` objects, including:
#' 
#' - **Abundance Data**: Retrieve raw, relative, or normalized abundance data with options
#'   for taxonomic rank-specific aggregation and output formatting.
#' - **Sample and Taxa Information**: Access metadata related to samples and taxa, such as sample
#'   IDs, taxa IDs, taxonomic classification, and metadata variables.
#' - **Network Data**: Extract network interaction data and community detection results, providing
#'   insights into the relationships and structure within microbial communities.
#'
#' Each getter function is designed to be intuitive and flexible, offering parameters that allow
#' for customized data retrieval based on the user's specific needs. Additionally, for collections
#' of `mgnet` objects managed within an `mgnetList`, these functions facilitate batch processing
#' across multiple datasets.
#'
#' @section Available Getter Functions:
#' - `abundance(object, rank, .fmt)`: Retrieves abundance data.
#' - `rel_abundance(object, rank, .fmt)`: Retrieves relative abundance data.
#' - `norm_abundance(object, rank, .fmt)`: Retrieves normalized abundance data.
#' - `sample_id(object)`: Retrieves sample IDs.
#' - `taxa_id(object)`: Retrieves taxa IDs.
#' - `info_sample(object, .fmt)`: Retrieves sample metadata.
#' - `info_taxa(object, .fmt)`: Retrieves taxa metadata.
#' - `network(object)`: Retrieves network interaction data.
#' - `community(object)`: Retrieves community detection results.
#'
#' @section Usage:
#' \dontrun{
#' # Retrieve raw abundance data as a data frame for a specified taxonomic rank
#' df_abundance <- abundance(mgnet_obj, rank = "Genus", .fmt = "df")
#'
#' # Get sample IDs from an mgnet object
#' sample_ids <- sample_id(mgnet_obj)
#'
#' # Access sample metadata as a tibble
#' sample_metadata <- info_sample(mgnet_obj, .fmt = "tbl")
#'
#' # Extract network data from an mgnet object
#' net_data <- network(mgnet_obj)
#'
#' # For mgnetList objects, apply a getter function across all contained mgnet objects
#' list_sample_metadata <- info_sample(mgnet_list, .fmt = "tbl")
#' }
#'
#' @name mgnet-getters
#' @rdname mgnet-getters
