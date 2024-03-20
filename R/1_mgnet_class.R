setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS MGNET
################################################################################
#' MetaGenomic NETwork (mgnet) Class
#'
#' An S4 class for comprehensive management and analysis of metagenomic networks.
#' It encapsulates various types of metagenomic data, including abundance matrices,
#' sample metadata, taxonomic information, and network analysis results.
#' This class serves as a robust framework for analyzing microbial communities
#' and their interactions, supporting advanced analysis workflows for metagenomics research.
#'
#' @slot abundance A numeric matrix representing next-generation sequencing (NGS)
#'        data, where rows correspond to samples and columns to taxa, with values indicating
#'        the abundance of each taxa in each sample.
#' @slot log_abundance A numeric matrix of log-ratio transformed abundance data,
#'        facilitating compositional data analysis, mirroring the structure of `abundance`.
#' @slot info_sample A data.frame containing metadata for samples, indexed by sample IDs,
#'        with each row corresponding to a sample and each column to an experimental variable.
#' @slot lineage A character matrix with taxonomic classification for each taxa,
#'        where each row corresponds to a taxa and columns to taxonomic ranks such as phylum,
#'        class, order, family, genus, and species.
#' @slot info_taxa A data.frame with additional information on taxa, indexed by taxa IDs,
#'        providing a flexible structure for storing diverse metadata related to taxa.
#' @slot network An `igraph` object representing a network of taxa interactions,
#'        capturing the complex relationships between microbial taxa within the community.
#' @slot community An object storing community detection results, typically obtained
#'        from network analysis, facilitating the exploration of microbial community structure.
#'
#' @section Reserved Keywords:
#' The `mgnet` class reserves two keywords for internal use, which should not be used
#' as column names in the provided matrices or data.frames: `sample_id`, `taxa_id`.
#' These keywords are utilized for specific functionalities and methods within the `mgnet` package.
#' Using these names as column identifiers in your data may lead to unexpected behavior or errors.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods setClass
#'
#' @seealso \code{\link{mgnet}} for the constructor function.
#'
#' @examples
#' # Creating an empty mgnet object
#' empty_mgnet <- mgnet()
#'
#' @name mgnet-class
#' @rdname mgnet-class
#' @exportClass mgnet
setClass("mgnet",
         
  slots = c(
    abundance = "ANY",
    log_abundance = "ANY",
    info_sample = "ANY",
    lineage = "ANY",
    info_taxa = "ANY",
    network = "ANY",
    community = "ANY"  
  ),
  
  prototype = prototype(
    abundance = matrix(numeric(0), nrow=0, ncol=0),
    log_abundance = matrix(numeric(0), nrow=0, ncol=0),
    info_sample = data.frame(),
    lineage = matrix(character(0), nrow=0, ncol=0),
    info_taxa = data.frame(),
    network = igraph::make_empty_graph(0),
    community = igraph::cluster_fast_greedy(
      igraph::make_empty_graph(0,directed=F))
  )
)
################################################################################
# END CLASS MGNET
################################################################################


################################################################################
# CONSTRUCTOR MGNET
################################################################################
#' Create an mgnet Object
#'
#' This function constructs an `mgnet` object, encapsulating various types of metagenomic data
#' including abundance data, sample metadata, taxonomic information, network analysis results, etc.,
#' providing a comprehensive framework for the analysis of microbial communities and their interactions.
#' It handles object creation with custom error handling for improved user experience.
#'
#' @param abundance numeric matrix with all elements >=0 representing the abundance data 
#'        where rows are samples and columns are taxa (OTUs, Species, ...). Defaults to an empty matrix.
#' @param log_abundance numeric matrix of log-transformed abundance data, structured like `abundance`.
#'        Defaults to an empty matrix.
#' @param info_sample data.frame containing metadata for samples. Each row should correspond
#'        to a sample and each column to an experimental variable. Defaults to an empty data frame.
#' @param lineage character matrix with taxonomic classification for each taxa. Each row
#'        corresponds to a taxa and each column to a rank. Defaults to an empty matrix.
#' @param info_taxa data.frame with additional information on taxa. Each row should correspond
#'        to a taxa. Defaults to an empty data frame.
#' @param network an `igraph` object representing a network of taxa interactions. Defaults to an
#'        empty graph.
#' @param community an object storing community detection results from `network`. Defaults
#'        to the result of a fast greedy clustering on an empty graph.
#'
#' @return Returns an `mgnet` object encapsulating the provided metagenomic data.
#'         If an error occurs during object creation, a custom error message is displayed and
#'         the function execution is stopped.
#'         
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods new
#'
#' @examples
#' # Assuming correct_matrix is a properly defined matrix and other necessary data are defined
#' mgnet_obj <- mgnet(abundance = otu_HMP2,
#'                    info_sample = info_sample_HMP2,
#'                    lineage = lineage_HMP2)
#' @export
mgnet <- function(abundance=matrix(nrow=0,ncol=0),
                  log_abundance=matrix(nrow=0, ncol=0),
                  info_sample=data.frame(),
                  lineage=matrix(nrow=0,ncol=0),
                  info_taxa=data.frame(),
                  network=make_empty_graph(n=0, directed=FALSE),
                  community=cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))
){
  
  if("sample_sum" %in% colnames(info_sample)) {
    stop(
      "The 'sample_sum' column should not be manually provided in 'info_sample'. ",
      "This column is essential for calculating relative abundances and is automatically calculated based on the provided 'abundance' data and it is generated during object construction.",
      "Alternatively, you can explicitly update 'sample_sum' using the 'update_sample_sum()' function, ",
    )
  }
  # Attempt to create a new mgnet object within a tryCatch block
  tryCatch({
    mgnet_object <- new("mgnet", 
                        abundance = abundance, log_abundance=log_abundance,
                        info_sample = info_sample,
                        lineage = lineage, info_taxa=info_taxa,
                        network = network, community = community
    )
    
    if(length(abundance)!=0 && length(info_sample)==0){
      info_sample <- data.frame("sample_sum"=rowSums(abundance))
    } else if( length(info_sample)!=0 & !("sample_sum"%in%colnames(info_sample))){
      info_sample$sample_sum <- rowSums(abundance)
    }
    mgnet_object@info_sample <- info_sample
    return(mgnet_object)
  }, error = function(e) {
    
    # Extract and split the error message into lines
    error_lines <- strsplit(e$message, "\n")[[1]]
    
    # Filter lines that start with "DEBUGGER" or "-"
    filtered_lines <- grep("^(DEBUGGER|-)", error_lines, value = TRUE)
    
    # Concatenate the filtered lines into a single string
    custom_error_msg <- paste(filtered_lines, collapse = "\n")
    
    # Throw a stop with the cleaned error message
    stop("mgnet class error message:\n", custom_error_msg, "\n")
  })
}
################################################################################
# END CONSTRUCTOR MGNET
###############################################################################
