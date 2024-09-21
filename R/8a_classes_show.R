# ################################################################################
# # SHOW MGNET
# ################################################################################
#' @importFrom utils head tail
setMethod("show", "mgnet", function(object) {
  cat("==== mgnet Object Summary ====\n")

  # General information
  cat("General Info:\n")
  cat(sprintf("  Samples: %d\n", nsample(object)))
  cat(sprintf("  Taxa: %d\n", ntaxa(object)))
  
  if( length(object@abun) > 0){
    zeroPercentage <- sum(object@abun == 0) / length(object@rela)
    cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
  } else if (length(object@rela) > 0) {
    zeroPercentage <- sum(object@rela == 0) / length(object@rela)
    cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
  } else {
    cat("  No abundance or rela available.\n")
  }
  

  # Sample metadata
  if( length(object@sample) > 0) {
    colNamesSample <- names(object@sample)
    sampleInfo <- paste0("  Sample Meta Info: ", toString(head(colNamesSample, 4)))
    if(length(colNamesSample) > 4) sampleInfo <- paste(sampleInfo, ", etc...")
    cat(sampleInfo, "\n")
  } else {
    cat("  No sample metadata available.\n")
  }

  # Taxa metadata
  if( length(object@taxa) > 0) {
    colNamesTaxa <- names(object@taxa)
    taxaMetaInfo <- paste0("  Taxa Meta Info: ", toString(head(colNamesTaxa, 4)))
    if(length(colNamesTaxa) > 4) taxaMetaInfo <- paste(taxaMetaInfo, ", etc...")
    cat(taxaMetaInfo, "\n")
  } else {
    cat("  No taxa metadata available.\n")
  }

  # Network information
  if( length(object@netw) > 0) {
    cat(sprintf("  Network: %d edges, ~%.2f%% density\n", igraph::ecount(object@netw), 
                                                          100*igraph::edge_density(object@netw) ))
  } else {
    cat("  No network data available.\n")
  }

  # Community information
  if( length(object@comm) > 0) {
    cat(sprintf("  Detected Communities: %d\n", max(as.numeric(names(sizes(object@comm)[sizes(object@comm) > 1]))) ))
    
    sizes <- igraph::sizes(object@comm)
    sizes <- sizes[sizes > 1]
    sizes_msg <- paste0("  Community Sizes: ", toString(head(tail(sizes, -1),6)))
    if(length(sizes) > 6) sizes_msg <- paste(sizes_msg, ", etc...")
    cat(sizes_msg, "\n")
    isolated_msg <- paste0("  Isolated Nodes: ", sum(igraph::sizes(object@comm)==1))
    cat(isolated_msg, "\n")
    
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
    
    if( length(mgnetObj@abun) > 0){
      zeroPercentage <- sum(mgnetObj@abun == 0) / length(mgnetObj@rela)
      cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
    } else if (length(mgnetObj@rela) > 0) {
      zeroPercentage <- sum(mgnetObj@rela == 0) / length(mgnetObj@rela)
      cat(sprintf("  Zeros Percentage: ~%.2f%%\n", 100 * zeroPercentage))
    } else {
      cat("  No abundance or rela available.\n")
    }
    
    # Display a simplified network and community info if available
    if (!is.null(mgnetObj@netw) && igraph::vcount(mgnetObj@netw) > 0) {
      cat(sprintf("  Network: %d edges, ~%.2f%% density\n", igraph::ecount(mgnetObj@netw), 
                                                            100*igraph::edge_density(mgnetObj@netw)))
    } else {
      cat("  No network data available.\n")
    }
    if (!is.null(mgnetObj@comm) && length(mgnetObj@comm) > 0) {
      cat(sprintf("  Detected Communities: %d (%.2f%% Isolated)\n", 
                  max(as.numeric(names(igraph::sizes(mgnetObj@comm)[igraph::sizes(mgnetObj@comm) > 1]))),
                  100*sum(igraph::sizes(mgnetObj@comm)==1) / igraph::vcount(mgnetObj@netw)))
    } else {
      cat("  No community detection results available.\n")
    }
  }
  
  cat("\n")
  cat("==== End of mgnetList Object Summary ====\n")
})