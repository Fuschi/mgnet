setOldClass("igraph")
setOldClass("communities")
################################################################################
# CLASS MGNET
################################################################################
#' Metagenomics Network Class
#'
#' An S4 class designed to manage and analyze metagenomics networks, including
#' abundance data, taxa lineage, sample and taxa metadata.
#'
#' @slot abundance numeric matrix representing the abundance data.
#' @slot sampleInfo data.frame containing metadata about samples.
#' @slot lineage matrix with the taxonomic lineage of each taxa.
#' @slot taxaInfo data.frame with metadata about taxa.
#' @slot association matrix representing relationships between taxa, such as correlations or conditional independencies.
#' @slot network object belong to 'igraph' class depicting the network of interactions among taxa.
#' @slot community A 'communities' object representing detected communities within the network.
#' @export 
#' @importFrom methods setClass setValidity
#' @importFrom igraph make_empty_graph cluster_fast_greedy
setClass("mgnet",
         
  slots = c(
    abundance = "ANY",
    sampleInfo = "ANY",
    lineage = "ANY",
    taxaInfo = "ANY",
    association = "ANY",
    network = "ANY",
    community = "ANY"  
  ),
  
  prototype = prototype(
    abundance = matrix(numeric(0), nrow=0, ncol=0),
    sampleInfo = data.frame(),
    lineage = matrix(character(0), nrow=0, ncol=0),
    taxaInfo = data.frame(),
    association = matrix(numeric(0), nrow=0, ncol=0),
    network = igraph::make_empty_graph(0),
    community = igraph::cluster_fast_greedy(
      igraph::make_empty_graph(0,directed=F))
  )
)



setValidity("mgnet", function(object) {
  errors <- character()
  
  # Check abundance slot
  if (length(object@abundance)>0){
    errorsAbund <- .assertAbundance(object@abundance,"abundance")
    errors <- c(errors,errorsAbund)
  }
  
  # Check sampleInfo slot
  if (length(object@sampleInfo)>0){
    errorsSampleInfo <- .assertMetaInfo(object@sampleInfo,"sampleInfo")
    errors <- c(errors,errorsSampleInfo)
  }
  
  # Check lineage slot
  if (length(object@lineage)>0){
    errorsLineage <- .assertLineage(object@lineage,"lineage")
    errors <- c(errors,errorsLineage)
  }
  
  # Check taxaInfo slot
  if (length(object@taxaInfo)>0){
    errorsTaxaInfo <- .assertMetaInfo(object@taxaInfo,"taxaInfo")
    errors <- c(errors,errorsTaxaInfo)
  }
  
  # Check association slot
  if (length(object@association)>0){
    errorsAssociation <- .assertAssociation(object@association,"association")
    errors <- c(errors,association)
  }
  
  # Aggregate and return all errors
  if (length(errors) > 0) {
    errors[1] <- paste("\n",errors[1],sep="")
    return(paste(errors, collapse = "\n"))
  }
  TRUE
})
################################################################################
# END CLASS MG
################################################################################


