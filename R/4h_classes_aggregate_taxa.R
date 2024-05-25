# AGGREGATE TAXA
#------------------------------------------------------------------------------#
#' Aggregate Data at Higher Taxonomic Levels
#' 
#' @description
#' This function aggregates abundance data and lineage at a higher taxonomic level specified 
#' by the user. It calculates the sum of data abundances over all taxa that map to the same 
#' higher-level group and removes any ambiguous levels from the taxonomy table.
#' 
#' @param object An `mgnet` object.
#' @param rank Taxonomic level chosen for aggregation.
#' 
#' @return An updated `mgnet` object with aggregated data and taxa at the specified taxonomic level.
#' The slots `info_taxa`, `network`, and `community` are lost after the aggregation process.
#' 
#' @export
#' @name aggregate_taxa
#' @aliases aggregate_taxa,mgnet,character-method aggregate_taxa,mgnetList,character-method
setGeneric("aggregate_taxa", function(object, rank) standardGeneric("aggregate_taxa"))

setMethod("aggregate_taxa", c("mgnet", "character"),
          function(object, rank) {
            if(length(object@lineage)==0) stop("lineage slot cannot be empty")
            if(!(rank %in% ranks(object))) stop("Rank must be one of these possible choices: ", toString(ranks(object)))
            
            data.aggregate <- mgnet::abundance(object, rank)
            
            # Aggregate taxa data
            taxa.aggregate <- object@lineage[, seq_len(match(rank, mgnet::ranks(object)))]
            taxa.aggregate <- taxa.aggregate[!duplicated(taxa.aggregate), ]
            rownames(taxa.aggregate) <- taxa.aggregate[, rank]
            taxa.aggregate <- taxa.aggregate[colnames(data.aggregate), ]
            
            # Display message if additional slots are present
            if (length(object@info_taxa) != 0 || length(object@network) != 0 || length(object@community) != 0) {
              message("The information regarding the slots 'info_taxa', 'network', and 'community' are lost after the function 'aggregate_taxa'.")
            }
            
            # Return updated mgnet object
            return(new("mgnet", abundance = as.matrix(data.aggregate), 
                       info_sample = object@info_sample, 
                       lineage = as.matrix(taxa.aggregate)))
          })

setMethod("aggregate_taxa", c("mgnetList", "character"),
          function(object, rank) {
            
            object@mgnets <- sapply(object@mgnets, function(x) aggregate_taxa(x, rank),
                                    simplify = FALSE, USE.NAMES = TRUE)
            return(object)
            
          })
#------------------------------------------------------------------------------#
