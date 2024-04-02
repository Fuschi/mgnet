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
  
  if( length(object@abundance) > 0){
    zeroPercentage <- sum(object@abundance == 0) / (nrow(object@abundance) * ncol(object@abundance))
    cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
  } else {
    cat("  No abundance available.\n")
  }
  

  # Sample metadata
  if( length(object@info_sample) > 0) {
    colNamesSample <- names(object@info_sample)
    sampleInfo <- paste0("  Sample Meta Info: ", toString(head(colNamesSample, 4)))
    if(length(colNamesSample) > 4) sampleInfo <- paste(sampleInfo, ", etc...")
    cat(sampleInfo, "\n")
  } else {
    cat("  No sample metadata available.\n")
  }

  # Taxonomic lineage
  if( length(object@lineage) > 0) {
    lineageInfo <- paste("  Taxonomic Ranks:", toString(colnames(object@lineage)))
    cat(lineageInfo, "\n")
  } else {
    cat("  No taxonomic lineage information available.\n")
  }

  # Taxa metadata
  if( length(object@info_taxa) > 0) {
    colNamesTaxa <- names(object@info_taxa)
    taxaMetaInfo <- paste0("  Taxa Meta Info: ", toString(head(colNamesTaxa, 4)))
    if(length(colNamesTaxa) > 4) taxaMetaInfo <- paste(taxaMetaInfo, ", etc...")
    cat(taxaMetaInfo, "\n")
  } else {
    cat("  No taxa metadata available.\n")
  }

  # Network information
  if( length(object@network) > 0) {
    cat(sprintf("  Network: %d edges, ~%.2f density\n", igraph::ecount(object@network), 
                                                        igraph::edge_density(object@network) ))
  } else {
    cat("  No network data available.\n")
  }

  # Community information
  if( length(object@community) > 0) {
    cat(sprintf("  Detected Communities: %d\n", max(igraph::membership(object@community))))
    
    sizes <- igraph::sizes(object@community)
    sizes_msg <- paste0("  Community Sizes:", toString(head(sizes, 6)))
    if(length(sizes) > 6) sizes_msg <- paste(sizes_msg, ", etc...")
    cat(sizes_msg, "\n")

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
  names <- names(object@mgnets)
  for (i in seq_along(object@mgnets)) {
    cat(sprintf("\n  -- mgnet Object `%s` --\n", names[i]))
    mgnetObj <- object@mgnets[[i]]
    cat(sprintf("  Samples: %d\n", nsample(mgnetObj)))
    cat(sprintf("  Taxa: %d\n", ntaxa(mgnetObj)))
    
    if(length(mgnetObj@abundance)!=0){
      zeroPercentage <- sum(mgnetObj@abundance == 0) / (nrow(mgnetObj@abundance) * ncol(mgnetObj@abundance))
      cat(sprintf("  Zeros Percentage: ~%.2f%%\n", zeroPercentage))
    } else {
      cat("  No abundance available.\n")
    }
    
    # Display a simplified network and community info if available
    if (!is.null(mgnetObj@network) && igraph::vcount(mgnetObj@network) > 0) {
      cat(sprintf("  Network: %d edges, ~%.2f density\n", igraph::ecount(mgnetObj@network), 
                                                          igraph::edge_density(mgnetObj@network)))
    } else {
      cat("  No network data available.\n")
    }
    if (!is.null(mgnetObj@community) && length(mgnetObj@community) > 0) {
      cat(sprintf("  Detected Communities: %d\n", max(igraph::membership(mgnetObj@community))))
    } else {
      cat("  No community detection results available.\n")
    }
  }
  
  cat("\n")
  cat("==== End of mgnetList Object Summary ====\n")
})