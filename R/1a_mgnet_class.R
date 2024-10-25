setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS MGNET
################################################################################
#' MetaGenomic NETwork (mgnet) Class
#'
#' An S4 class designed for comprehensive management and analysis of metagenomic networks.
#' It encapsulates a variety of microbial data including abundance matrices, sample metadata,
#' taxonomic information, and results from network analysis. This class serves as a robust
#' framework for exploring microbial communities and their interactions, supporting
#' sophisticated workflows in metagenomics research. Designed to be compatible with the tidyverse,
#' `mgnet` facilitates seamless integration with pipelines that leverage tidyverse packages
#' for data manipulation, analysis, and visualization.
#'
#' @slot abun A numeric matrix representing the raw abundance data from
#'        next-generation sequencing (NGS), where rows correspond to samples and columns
#'        to taxa.
#' @slot rela A numeric matrix representing the relative abundance of each taxon
#'        within each sample. It is calculated as the proportion of each taxon's abundance
#'        relative to the total abundance in the sample, providing insights into the
#'        relative presence of taxa across samples.
#' @slot norm A numeric matrix of normalized abundance data, suitable for comparative
#'        analysis across samples. This includes transformations such as log-ratio transformation,
#'        rarefaction, or DESeq normalization, to address issues like compositional data biases.
#' @slot meta A data frame containing metadata for samples, indexed by sample IDs.
#'        Metadata can include sample collection details, experimental conditions,
#'        and clinical outcomes.
#' @slot taxa A data frame containing a broad range of metadata associated with each taxon,
#'        indexed by taxa IDs. This metadata may encompass any relevant ecological, genetic,
#'        or phenotypic information pertinent to each taxon.
#' @slot netw An `igraph` object that depicts the network of interactions among taxa
#'        based on their co-occurrence and other ecological metrics. This network can
#'        be used to infer ecological relationships and community structure.
#' @slot comm An object storing results from community detection analysis performed
#'        on the network. This can help identify clusters or groups of taxa that frequently
#'        interact or share similar ecological niches.
#'
#' @section Reserved Keywords:
#' The `mgnet` class reserves the following keywords for internal use, which should not be used
#' as column names in provided matrices or data frames: `sample_id`, `taxa_id`, `comm_id`,
#' `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`, `.`, `mgnet`.
#' These keywords are integral to the functioning of methods within the `mgnet` package, and using these as column identifiers cause errors.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods setClass
#'
#' @seealso \link{mgnet} for the constructor function which details how to create an instance
#'          of this class.
#'
#' @examples
#' # Creating an empty mgnet object
#' empty_mgnet <- mgnet()
#'
#' # Creating a mgnet object with example data
#' data(otu_HMP2, meta_HMP2, taxa_HMP2, package = "mgnet")
#' HMP2 <- mgnet(abun = otu_HMP2, meta = meta_HMP2, taxa = taxa_HMP2)
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
    abun = matrix(numeric(0), nrow = 0, ncol = 0),
    rela = matrix(numeric(0), nrow = 0, ncol = 0),
    norm = matrix(numeric(0), nrow = 0, ncol = 0),
    meta = data.frame(),
    taxa = data.frame(),
    netw = igraph::make_empty_graph(0),
    comm = igraph::cluster_fast_greedy(
      igraph::make_empty_graph(0, directed = F)
    )
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
#' This function constructs an `mgnet` object, encapsulating a comprehensive range of metagenomic data.
#' It integrates abundance data, sample metadata, taxonomic details, and network analysis results into a single
#' structure, providing a robust platform for analyzing microbial communities and their interactions.
#' Enhanced error handling is implemented to ensure a smooth user experience by providing clear feedback
#' on input errors.
#'
#' @details
#' If an abundance matrix, `abun`, is provided and the relative abundance matrix is omitted, the constructor
#' automatically calculates relative abundances and store them in `rela`.
#'
#' @param abun Numeric matrix representing the abundance of taxa across samples, where rows correspond to
#'        samples and columns to taxa (e.g., OTUs, species). Each element must be >=0, representing the
#'        abundance count of a taxon in a sample. Defaults to an empty matrix.
#' @param rela Numeric matrix representing the relative abundance of each taxon within each sample,
#'        calculated as the proportion of a taxon's abundance relative to the total abundance in the sample.
#'        Defaults to an empty matrix.
#' @param norm Numeric matrix of normalized abundance data, suitable for comparative analyses across samples,
#'        implementing normalization methods such as log-ratio transformations. Defaults to an empty matrix.
#' @param meta Data frame containing metadata for samples. Each row should correspond to a sample,
#'        with columns representing various experimental variables. Defaults to an empty data frame.
#' @param taxa Data frame containing a broad range of metadata associated with each taxon, indexed by taxa IDs,
#'        facilitating flexible integration of diverse ecological and genetic information. Defaults to an empty data frame.
#' @param netw An `igraph` object representing the network of interactions among taxa, which can be used to
#'        infer ecological relationships and community structures. Defaults to an empty graph created with
#'        `make_empty_graph`.
#' @param comm An object storing community detection results, typically obtained from network analysis.
#'        Defaults to the result of `cluster_fast_greedy` applied to an empty graph.
#'
#' @section Reserved Keywords:
#' The `mgnet` class reserves the following keywords for internal use, which should not be used
#' as column names in provided matrices or data frames: `sample_id`, `taxa_id`, `comm_id`,
#' `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`, `.`, `mgnet`.
#' These keywords are integral to the functioning of methods within the `mgnet` package, and using these as column identifiers cause errors.
#'
#' @return An `mgnet` object that encapsulates the provided metagenomic data. Errors during object creation
#'         trigger custom error messages and halt function execution, ensuring data integrity and user guidance.
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods new
#'
#' @examples
#' # Creating an empty mgnet object
#' empty_mgnet <- mgnet()
#'
#' # Creating a mgnet object with example data
#' data(otu_HMP2, meta_HMP2, taxa_HMP2, package = "mgnet")
#' HMP2 <- mgnet(abun = otu_HMP2, meta = meta_HMP2, taxa = taxa_HMP2)
#'
#' @export
mgnet <- function(abun = matrix(numeric(0), nrow = 0, ncol = 0),
                  rela = matrix(numeric(0), nrow = 0, ncol = 0),
                  norm = matrix(numeric(0), nrow = 0, ncol = 0),
                  meta = data.frame(),
                  taxa = data.frame(),
                  netw = make_empty_graph(n = 0, directed = FALSE),
                  comm = cluster_fast_greedy(make_empty_graph(n = 0, directed = FALSE))) {
  tryCatch(
    {
      mgnet_object <- new("mgnet",
        abun = abun,
        rela = rela, norm = norm,
        meta = meta,
        taxa = taxa,
        netw = netw, comm = comm
      )

      if (length(mgnet_object@abun) != 0 && length(mgnet_object@rela) == 0) {
        mgnet_object@rela <- mgnet_object@abun / rowSums(mgnet_object@abun)
      }

      return(mgnet_object)
    },
    error = function(e) {
      # Clean the error message to remove R internal trace information
      stop(e$message)
    }
  )
}
################################################################################
# END CONSTRUCTOR MGNET
################################################################################
