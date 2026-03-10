#' @include class-mgnet.R class-grouping.R class-links.R
NULL

# mgnet EXTRACTOR
#------------------------------------------------------------------------------#
#' Subset mgnet Object
#'
#' Subsets an `mgnet` object based on sample (i) and taxa (j) indices.
#'
#' @name subset-mgnet
#' @rdname subset-mgnet
#' @aliases [,mgnet-method \S4method{[}{mgnet}
#'
#' @usage \S4method{[}{mgnet}(x, i, j, ..., drop = FALSE)
#'
#' @param x An `mgnet` object to be subsetted.
#' @param i Indices or a logical vector selecting **samples**. If missing, all samples are included.
#' @param j Indices or a logical vector selecting **taxa**. If missing, all taxa are included.
#' @param ... Ignored. Present only for compatibility with the S4 `[` signature.
#' @param drop Ignored. Always treated as `FALSE`; passing `drop = TRUE` has no effect.
#'
#' @details
#' Subsetting by samples filters abundance matrices and sample metadata; subsetting by taxa
#' filters taxa-related data and, if present, updates dependent structures accordingly.
#' Subsetting by samples removes the `netw` and `comm` slots; a warning is issued.
#'
#' @return A new `mgnet` object containing only the requested subset.
#' @seealso \code{\link[base]{Extract}}
#' @importFrom igraph make_empty_graph cluster_fast_greedy
#' @importFrom methods new
#' @importFrom igraph V E induced_subgraph permute membership make_clusters
#' @importFrom cli cli_abort cli_warn
#' @export
setMethod(f = "[", signature = "mgnet", definition = function(x, i, j, ..., drop = FALSE) {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  check_unique <- function(indices, label) {
    if (any(duplicated(indices))) {
      cli::cli_abort("'{.var {label}}' contains duplicate values.")
    }
  }
  
  empty_comm <- function() {
    igraph::cluster_fast_greedy(
      igraph::make_empty_graph(0, directed = FALSE)
    )
  }
  
  # Validate and process 'i' indices (samples)
  if (!missing(i)) {
    if (is.numeric(i)) {
      if (any(i < 1 | i > nsample(x))) {
        cli::cli_abort("{.var i} (samples) indices out of range (1..{nsample(x)}).")
      }
      check_unique(i, "i (samples)")
    } else if (is.character(i)) {
      i <- match(i, sample_id(x))
      if (any(is.na(i))) {
        cli::cli_abort("Some sample names in {.var i} not found in {.fn sample_id}().")
      }
      check_unique(i, "i (samples)")
    } else if (is.logical(i)) {
      if (length(i) != nsample(x)) {
        cli::cli_abort("Logical {.var i} must have length equal to the number of samples ({nsample(x)}).")
      }
      i <- which(i)
    } else {
      cli::cli_abort("Invalid {.var i} (samples). Must be numeric, character, or logical.")
    }
  } else {
    i <- seq_len(nsample(x))
  }
  
  # Validate and process 'j' indices (taxa)
  if (!missing(j)) {
    if (is.numeric(j)) {
      if (any(j < 1 | j > ntaxa(x))) {
        cli::cli_abort("{.var j} (taxa) indices out of range (1..{ntaxa(x)}).")
      }
      check_unique(j, "j (taxa)")
    } else if (is.character(j)) {
      j <- match(j, taxa_id(x))
      if (any(is.na(j))) {
        cli::cli_abort("Some taxa names in {.var j} not found in {.fn taxa_id}().")
      }
      check_unique(j, "j (taxa)")
    } else if (is.logical(j)) {
      if (length(j) != ntaxa(x)) {
        cli::cli_abort("Logical {.var j} must have length equal to the number of taxa ({ntaxa(x)}).")
      }
      j <- which(j)
    } else {
      cli::cli_abort("Invalid {.var j} (taxa). Must be numeric, character, or logical.")
    }
  } else {
    j <- seq_len(ntaxa(x))
  }
  
  # Full sample set retained?
  full_i <- identical(sort(i), seq_len(nsample(x)))
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # Handle empty subset cases
  if (length(i) == 0 && length(j) == 0) {
    
    new_mgnet <- new("mgnet")
    
  } else if (length(i) == 0 && length(j) > 0) {
    
    # Case where only taxa are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(numeric(0), nrow = 0, ncol = 0)
    new_mgnet@rela <- matrix(numeric(0), nrow = 0, ncol = 0)
    new_mgnet@norm <- matrix(numeric(0), nrow = 0, ncol = 0)
    new_mgnet@meta <- data.frame()
    
    if (length(x@taxa) != 0) {
      new_mgnet@taxa <- x@taxa[j, , drop = FALSE]
    } else {
      new_mgnet@taxa <- x@taxa
    }
    
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm <- empty_comm()
    
  } else if (length(i) > 0 && length(j) == 0) {
    
    # Case where only samples are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(numeric(0), nrow = 0, ncol = 0)
    new_mgnet@rela <- matrix(numeric(0), nrow = 0, ncol = 0)
    new_mgnet@norm <- matrix(numeric(0), nrow = 0, ncol = 0)
    
    if (length(x@meta) != 0) {
      new_mgnet@meta <- x@meta[i, , drop = FALSE]
    } else {
      new_mgnet@meta <- x@meta
    }
    
    new_mgnet@taxa <- data.frame()
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm <- empty_comm()
    
  } else {
    
    # Standard case where both i and j are present
    if (length(x@abun) != 0) abun.new <- x@abun[i, j, drop = FALSE] else abun.new <- x@abun
    if (length(x@rela) != 0) rela.new <- x@rela[i, j, drop = FALSE] else rela.new <- x@rela
    if (length(x@norm) != 0) norm.new <- x@norm[i, j, drop = FALSE] else norm.new <- x@norm
    if (length(x@meta) != 0) meta.new <- x@meta[i, , drop = FALSE] else meta.new <- x@meta
    if (length(x@taxa) != 0) taxa.new <- x@taxa[j, , drop = FALSE] else taxa.new <- x@taxa
    
    if (length(x@netw) == 0) {
      
      # No network: just subset abundances and metadata
      new_mgnet <- mgnet(
        abun = abun.new,
        rela = rela.new,
        norm = norm.new,
        meta = meta.new,
        taxa = taxa.new
      )
      
    } else if (length(x@netw) != 0 && full_i) {
      
      # Network can be retained only if all samples are preserved
      netw <- x@netw
      netw.new <- igraph::induced_subgraph(netw, vids = igraph::V(netw)[j])
      netw.new <- igraph::permute(
        netw.new,
        match(igraph::V(netw.new)$name, igraph::V(netw)$name[j])
      )
      
      if (length(x@comm) != 0) {
        v_old <- igraph::V(netw)$name
        v_new <- igraph::V(netw.new)$name
        
        m_new <- igraph::membership(x@comm)[match(v_new, v_old)]
        names(m_new) <- v_new
        
        algo <- attr(x@comm, "algorithm")
        if (is.null(algo)) algo <- "subset"
        
        comm.new <- igraph::make_clusters(
          netw.new,
          membership = m_new,
          algorithm  = algo,
          modularity = NA_real_
        )
      } else {
        comm.new <- x@comm
      }
      
      new_mgnet <- mgnet(
        abun = abun.new,
        rela = rela.new,
        norm = norm.new,
        meta = meta.new,
        taxa = taxa.new,
        netw = netw.new,
        comm = comm.new
      )
      
      # Preserve selected links by surviving link_id values
      if (are_selected_links(x)) {
        old_selection <- get_selected_links(x)
        new_selection <- old_selection[old_selection %in% igraph::E(netw.new)$link_id]
        attr(new_mgnet, "selected_links") <- new_selection
      }
      
      # Preserve edge-level link grouping by surviving link_id values
      if (is_link_grouped(x)) {
        old_groups <- get_grouped_link(x)
        old_ids    <- igraph::E(x@netw)$link_id
        new_ids    <- igraph::E(netw.new)$link_id
        
        idx <- match(new_ids, old_ids)
        new_grouping <- old_groups[idx]
        
        attr(new_mgnet, "link_groups") <- new_grouping
      }
      
    } else if (length(x@netw) != 0 && !full_i) {
      
      # Subsetting by samples invalidates network/community coherence
      cli::cli_warn("Subsetting by samples removes {.field netw} and {.field comm} slots.")
      
      new_mgnet <- mgnet(
        abun = abun.new,
        rela = rela.new,
        norm = norm.new,
        meta = meta.new,
        taxa = taxa.new
      )
      
    } else {
      
      cli::cli_abort("Unexpected subsetting state in [.mgnet].")
    }
  }
  
  # Reapply meta/taxa grouping if it exists and is still representable
  if (is_mgnet_grouped(x)) {
    new_mgnet <- group_mgnet(new_mgnet, !!!rlang::syms(get_group_mgnet(x)))
  }
  
  new_mgnet
})