setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS mgnet
################################################################################
#' Class mgnet
#'
#' An S4 class designed for organizing and managing NGS data in a tidy and modular fashion.
#' The class stores key components such as abundance matrices (raw, relative, normalized),
#' sample metadata, feature annotations, network structures, and community assignments.
#'
#' This class is designed to integrate NGS data — abundances, metadata, annotations, networks —
#' with structural consistency checks. It performs no modeling or normalization and focuses on
#' tidy-style manipulation and coherent data export.
#'
#' @slot abun A numeric matrix of raw abundance values (samples x features).
#' @slot rela A numeric matrix of relative abundances.
#' @slot norm A numeric matrix of normalized abundances.
#' @slot meta A data.frame with sample metadata.
#' @slot taxa A data.frame with feature (e.g. taxon, gene) metadata.
#' @slot netw An igraph object representing a network among features or samples.
#' @slot comm A communities object (from igraph) representing community structure.
#'
#' @section Reserved Keywords:
#' The following reserved names are used internally and must be avoided:
#' `r paste(.MGNET_RESERVED_COLNAMES, collapse = ", ")`
#'
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods setClass
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
           netw = igraph::make_empty_graph(0, directed = F),
           comm = igraph::cluster_fast_greedy(
             igraph::make_empty_graph(0, directed = F)
           )
         )
)


################################################################################
# CONSTRUCTOR mgnet
################################################################################
#' Create an mgnet Object
#'
#' This constructor creates an `mgnet` object designed to encapsulate a single
#' microbial community dataset derived from high-throughput sequencing techniques.
#' The object integrates raw, relative and normalized abundance data, sample and taxa metadata, and
#' a taxon interaction network and its community structure.
#'
#' @param abun A numeric matrix of raw abundances (e.g., read counts), where rows
#'   represent samples and columns represent features (taxa, genes, etc.).
#'   All values must be ≥ 0. Default is an empty matrix.
#' @param rela A numeric matrix of relative abundances, where each value represents
#'   the proportion of a feature in a given sample. Defaults to empty. 
#'   All values must be ≥ 0. Default is an empty matrix.
#' @param norm A numeric matrix of transformed or normalized abundances, such as
#'   log(CPM), CLR-transformed data, or TPM. Defaults to an empty matrix.
#' @param meta A data frame of sample metadata. Defaults an empty data.frame.
#' @param taxa A data frame of taxon metadata. Defaults an empty data.frame.
#' @param netw An `igraph` object representing a feature-level network. Defaults to an empty graph.
#' @param comm A `communities` object representing detected communities in the network.
#'   Defaults to an empty result based on `cluster_fast_greedy`.
#'
#' @return An `mgnet` object. If invalid, the constructor halts with an informative error.
#'
#' @section Reserved Keywords:
#' The following reserved names are used internally and must be avoided:
#' `r paste(.MGNET_RESERVED_COLNAMES, collapse = ", ")`
#'
#' @examples
#' # Construct an mgnet object with example data from HMP2
#' HMP2 <- mgnet(
#'   abun = otu_HMP2,
#'   meta = meta_HMP2,
#'   taxa = taxa_HMP2
#' )
#'
#' @importFrom methods new
#' @export
#' @rdname mgnet-class
mgnet <- function(abun = matrix(numeric(0), nrow = 0, ncol = 0),
                 rela = matrix(numeric(0), nrow = 0, ncol = 0),
                 norm = matrix(numeric(0), nrow = 0, ncol = 0),
                 meta = data.frame(),
                 taxa = data.frame(),
                 netw = igraph::make_empty_graph(n = 0, directed = FALSE),
                 comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(n = 0, directed = FALSE))) {
  
  # --- Ensure link_id exists if edges are present ----------------------------#
  if(length(netw) > 0) netw <- .ensure_link_id(netw)
  
  # --- Construct the object --------------------------------------------------#
  tryCatch({
    methods::new("mgnet",
        abun = abun,
        rela = rela,
        norm = norm,
        meta = meta,
        taxa = taxa,
        netw = netw,
        comm = comm
    )
  }, error = function(e) {
    cli::cli_abort("Failed to create {.cls mgnet} object:\n{e$message}")
  })
}
