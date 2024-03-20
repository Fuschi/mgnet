# ################################################################################
# # SHOW MGNET
# ################################################################################
#' @importFrom utils head
setMethod("show", "mgnet", function(object) {
  cat("==== mgnet Object Summary ====\n")

  # General information
  cat("General Info:\n")
  cat(sprintf("  Samples: %d\n", nsample(object)))
  cat(sprintf("  Taxa: %d\n", ntaxa(object)))
  
  if(length(object@abundance)!=0){
    zeroPercentage <- sum(object@abundance == 0) / (nrow(object@abundance) * ncol(object@abundance))
    cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
  } else {
    cat("  No abundance available.\n")
  }
  

  # Sample metadata
  if(!is.null(object@info_sample) && ncol(object@info_sample) > 0) {
    colNamesSample <- names(object@info_sample)
    sampleInfo <- paste0("  Sample Meta Info: ", toString(head(colNamesSample, 4)))
    if(length(colNamesSample) > 4) sampleInfo <- paste(sampleInfo, ", etc...")
    cat(sampleInfo, "\n")
  } else {
    cat("  No sample metadata available.\n")
  }

  # Taxonomic lineage
  if(!is.null(object@lineage) && ncol(object@lineage) > 0) {
    lineageInfo <- paste("  Taxonomic Ranks:", toString(colnames(object@lineage)))
    cat(lineageInfo, "\n")
  } else {
    cat("  No taxonomic lineage information available.\n")
  }

  # Taxa metadata
  if(!is.null(object@info_taxa) && ncol(object@info_taxa) > 0) {
    colNamesTaxa <- names(object@info_taxa)
    taxaMetaInfo <- paste0("  Taxa Meta Info: ", toString(head(colNamesTaxa, 4)))
    if(length(colNamesTaxa) > 4) taxaMetaInfo <- paste(taxaMetaInfo, ", etc...")
    cat(taxaMetaInfo, "\n")
  } else {
    cat("  No taxa metadata available.\n")
  }

  # Network information
  if(!is.null(object@network) && igraph::vcount(object@network) > 0) {
    cat(sprintf("  Network: %d nodes, %d edges\n", igraph::vcount(object@network), igraph::ecount(object@network)))
    density <- igraph::edge_density(object@network)
    cat(sprintf("  Edge Density: %.4f\n", density))
  } else {
    cat("  No network data available.\n")
  }

  # Community information
  if(!is.null(object@community) && length(object@community) > 0) {
    cat(sprintf("  Detected Communities: %d\n", max(igraph::membership(object@community))))
    sizes <- toString(igraph::sizes(object@community))
    cat(sprintf("  Community Sizes: %s\n", sizes))
  } else {
    cat("  No community detection results available.\n")
  }

  cat("==== End of mgnet Object Summary ====\n")
})


# ################################################################################
# # SHOW MGNETLIST
# ################################################################################
#' @importFrom utils head
setMethod("show", "mgnetList", function(object) {
  cat("==== mgnetList Object Summary ====\n")
  
  # Number of mgnet objects in the list
  cat(sprintf("Contains %d mgnet objects:\n", length(object@mgnets)))
  
  # Iterate through each mgnet object and provide a summary
  for (i in seq_along(object@mgnets)) {
    cat(sprintf("\n  -- mgnet Object %d --\n", i))
    mgnetObj <- object@mgnets[[i]]
    cat(sprintf("  Samples: %d\n", nsample(mgnetObj)))
    cat(sprintf("  Taxa: %d\n", ntaxa(mgnetObj)))
    
    if(length(mgnetObj@abundance)!=0){
      zeroPercentage <- sum(mgnetObj@abundance == 0) / (nrow(mgnetObj@abundance) * ncol(mgnetObj@abundance))
      cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
    } else {
      cat("  No abundance available.\n")
    }
    
    # Display a simplified network and community info if available
    if (!is.null(mgnetObj@network) && igraph::vcount(mgnetObj@network) > 0) {
      cat(sprintf("  Network: %d nodes, %d edges\n", igraph::vcount(mgnetObj@network), igraph::ecount(mgnetObj@network)))
    } else {
      cat("  No network data available.\n")
    }
    if (!is.null(mgnetObj@community) && length(mgnetObj@community) > 0) {
      cat(sprintf("  Detected Communities: %d\n", max(igraph::membership(mgnetObj@community))))
    } else {
      cat("  No community detection results available.\n")
    }
  }
  
  cat("==== End of mgnetList Object Summary ====\n")
})