setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS MGNET
################################################################################
#' MetaGenomic NETwork (mgnet) Class
#'
#' An S4 class for comprehensive management and analysis of metagenomic networks.
#' It encapsulates various types of microbial data, including abundance matrices,
#' sample metadata, taxonomic information, and network analysis results.
#' This class serves as a robust framework for analyzing microbial communities
#' and their interactions, supporting advanced analysis workflows for metagenomics research.
#'
#' @slot abundance A numeric matrix representing the raw abundance data from
#'        next-generation sequencing (NGS), where rows correspond to samples and columns
#'        to taxa. Each entry indicates the abundance of a taxon in a sample.
#' @slot rel_abundance A numeric matrix representing the relative abundance of each taxon
#'        within each sample, mirroring the structure of `abundance`. This is calculated as
#'        the proportion of each taxon's abundance relative to the total abundance in the sample.
#' @slot norm_abundance A numeric matrix of normalized abundance data, facilitating various
#'        types of analyses. While the default normalization applied is log-ratio transformation
#'        to address compositional data issues, this slot is versatile and can hold other types
#'        of normalized data, such as those obtained from rarefaction, DESeq normalization, etc.
#' @slot info_sample A data.frame containing metadata for samples, indexed by sample IDs.
#' @slot lineage A character matrix providing taxonomic classification for each taxon.
#' @slot info_taxa A data.frame with additional taxa information, indexed by taxa IDs.
#' @slot network An `igraph` object representing the network of taxa interactions.
#' @slot community An object storing community detection results from network analysis.
#'
#'
#' @section Reserved Keywords:
#' The `mgnet` class reserves keywords for internal use, which should not be used
#' as column names in the provided matrices or data.frames: `sample_id`, `taxa_id`,
#' `abundance`,`rel_abundance`, `norm_abundance`.
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
    rel_abundance = "ANY",
    norm_abundance = "ANY",
    info_sample = "ANY",
    lineage = "ANY",
    info_taxa = "ANY",
    network = "ANY",
    community = "ANY"  
  ),
  
  prototype = prototype(
    abundance = matrix(numeric(0), nrow=0, ncol=0),
    rel_abundance = matrix(numeric(0), nrow=0, ncol=0),
    norm_abundance = matrix(numeric(0), nrow=0, ncol=0),
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
#' @param abundance Numeric matrix with all elements >=0 representing the abundance data 
#'        where rows are samples and columns are taxa (OTUs, Species, ...). Defaults to an empty matrix.
#' @param rel_abundance Numeric matrix representing the relative abundance of each taxon within each sample.
#'        Defaults to an empty matrix.
#' @param norm_abundance Numeric matrix of log-transformed abundance data, structured like `abundance`.
#'        Defaults to an empty matrix.
#' @param info_sample Data.frame containing metadata for samples. Each row should correspond
#'        to a sample and each column to an experimental variable. Defaults to an empty data frame.
#' @param lineage Character matrix with taxonomic classification for each taxa. Each row
#'        corresponds to a taxa and each column to a rank. Defaults to an empty matrix.
#' @param info_taxa Data.frame with additional information on taxa. Each row should correspond
#'        to a taxa. Defaults to an empty data frame.
#' @param network An `igraph` object representing a network of taxa interactions. Defaults to an
#'        empty graph.
#' @param community An object storing community detection results from `network`. Defaults
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
mgnet <- function(abundance = matrix(numeric(0), nrow=0,ncol=0),
                  rel_abundance = matrix(numeric(0), nrow=0,ncol=0),
                  norm_abundance = matrix(numeric(0), nrow=0, ncol=0),
                  info_sample = data.frame(),
                  lineage = matrix(character(0), nrow=0,ncol=0),
                  info_taxa = data.frame(),
                  network = make_empty_graph(n=0, directed=FALSE),
                  community = cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))
){
  
  tryCatch({
    mgnet_object <- new("mgnet", 
                        abundance = abundance, 
                        rel_abundance = rel_abundance, norm_abundance = norm_abundance,
                        info_sample = info_sample,
                        lineage = lineage, info_taxa = info_taxa,
                        network = network, community = community)
    
    if(length(mgnet_object@abundance) != 0 && length(mgnet_object@rel_abundance)==0 ){
      mgnet_object@rel_abundance <- mgnet_object@abundance / rowSums(mgnet_object@abundance)
    }
    
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
################################################################################
