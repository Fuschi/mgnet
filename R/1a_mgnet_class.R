setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS MGNET
################################################################################
#' MetaGenomic NETwork (mgnet) Class
#'
#' An S4 class for comprehensive management and analysis of metagenomic networks.
#' It encapsulates various types of microbial data, including abundances matrices,
#' sample metadata, taxonomic information, and network analysis results.
#' This class serves as a robust framework for analyzing microbial communities
#' and their interactions, supporting advanced analysis workflows for metagenomics research.
#'
#' @slot abun A numeric matrix representing the raw abundance data from
#'        next-generation sequencing (NGS), where rows correspond to samples and columns
#'        to taxa. Each entry indicates the abundance of a taxon in a sample.
#' @slot rela A numeric matrix representing the relative abundance of each taxon
#'        within each sample, mirroring the structure of `abundance`. This is calculated as
#'        the proportion of each taxon's abundance relative to the total abundance in the sample.
#' @slot norm A numeric matrix of normalized abundance data, facilitating various
#'        types of analyses. While the default normalization applied is log-ratio transformation
#'        to address compositional data issues, this slot is versatile and can hold other types
#'        of normalized data, such as those obtained from rarefaction, DESeq normalization, etc.
#' @slot meta A data.frame containing metadata for samples, indexed by sample IDs.
#' @slot taxa A data.frame with additional taxa information, indexed by taxa IDs.
#' @slot netw An `igraph` object representing the network of taxa interactions.
#' @slot comm An object storing community detection results from network analysis.
#'
#'
#' @section Reserved Keywords:
#' The `mgnet` class reserves keywords for internal use, which should not be used
#' as column names in the provided matrices or data.frames: `sample_id`, `taxa_id`, `comm_id`,
#' `abun`,`rela`, `norm`.
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
    abun = "ANY",
    rela = "ANY",
    norm = "ANY",
    meta = "ANY",
    taxa = "ANY",
    netw = "ANY",
    comm = "ANY"  
  ),
  
  prototype = prototype(
    abun = matrix(numeric(0), nrow=0, ncol=0),
    rela = matrix(numeric(0), nrow=0, ncol=0),
    norm = matrix(numeric(0), nrow=0, ncol=0),
    meta = data.frame(),
    taxa = data.frame(),
    netw = igraph::make_empty_graph(0),
    comm = igraph::cluster_fast_greedy(
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
#' including abundances data, sample metadata, taxonomic information, network analysis results, etc.,
#' providing a comprehensive framework for the analysis of microbial communities and their interactions.
#' It handles object creation with custom error handling for improved user experience.
#'
#' @param abun Numeric matrix with all elements >=0 representing the abundance data 
#'        where rows are samples and columns are taxa (OTUs, Species, ...). Defaults to an empty matrix.
#' @param rela Numeric matrix representing the relative abundance of each taxon within each sample.
#'        Defaults to an empty matrix.
#' @param norm Numeric matrix of normalized abundance data.
#'        Defaults to an empty matrix.
#' @param meta Data.frame containing metadata for samples. Each row should correspond
#'        to a sample and each column to an experimental variable. Defaults to an empty data frame.
#' @param taxa Data.frame with additional information on taxa. Each row should correspond
#'        to a taxa. Defaults to an empty data frame.
#' @param netw An `igraph` object representing a network of taxa interactions. Defaults to an
#'        empty graph.
#' @param comm An object storing community detection results from `network`. Defaults
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
#' mgnet_obj <- mgnet(abun = otu_HMP2,
#'                    sample = info_sample_HMP2,
#'                    taxa = info_taxa_HMP2)
#' @export
mgnet <- function(abun = matrix(numeric(0), nrow=0,ncol=0),
                  rela = matrix(numeric(0), nrow=0,ncol=0),
                  norm = matrix(numeric(0), nrow=0, ncol=0),
                  meta = data.frame(),
                  taxa = data.frame(),
                  netw = make_empty_graph(n=0, directed=FALSE),
                  comm = cluster_fast_greedy(make_empty_graph(n=0, directed=FALSE))
){
  
  tryCatch({
    mgnet_object <- new("mgnet", 
                        abun = abun, 
                        rela = rela, norm = norm,
                        meta = meta,
                        taxa = taxa,
                        netw = netw, comm = comm)
    
    if(length(mgnet_object@abun) != 0 && length(mgnet_object@rela)==0 ){
      mgnet_object@rela <- mgnet_object@abun / rowSums(mgnet_object@abun)
    }
    
    return(mgnet_object)
    
  }, error = function(e) {
    # Clean the error message to remove R internal trace information
    stop(e$message)
  })
}
################################################################################
# END CONSTRUCTOR MGNET
################################################################################
