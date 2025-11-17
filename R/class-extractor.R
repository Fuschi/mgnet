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
#' Subsetting by samples removes the `network` and `community` slots; a warning is issued.
#'
#' @return A new `mgnet` object containing only the requested subset.
#' @seealso \code{\link[base]{Extract}}
#' @importFrom igraph make_empty_graph cluster_fast_greedy 
#' @importFrom methods new 
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix V induced_subgraph
#' @importFrom cli cli_abort cli_warn
#' @export
setMethod(f="[", signature="mgnet", definition = function(x, i, j, ..., drop = FALSE)  {
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Helper function to check uniqueness (only for numeric or character indices)
  check_unique <- function(indices, label) {
    if (any(duplicated(indices))) {
      cli::cli_abort("'{.var {label}}' contains duplicate values.")
    }
  }
  
  # Validate and process 'i' indices (samples)
  if (!missing(i)) {
    if (is.numeric(i)) {
      if (any(i < 1 | i > nsample(x))) {
        cli::cli_abort("{.var i} (samples) indices out of range (1..{nsample(x)}).")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.character(i)) {
      i <- match(i, sample_id(x))
      if (any(is.na(i))) {
        cli::cli_abort("Some sample names in {.var i} not found in {.fn sample_id}().")
      }
      check_unique(i, "i (samples)")  # Check uniqueness
    } else if (is.logical(i)) {
      if (length(i) != nsample(x)) {
        cli::cli_abort("Logical {.var i} must have length equal to the number of samples ({nsample(x)}).")
      }
      i <- which(i)  # Convert logical to numeric
    } else {
      cli::cli_abort("Invalid {.var i} (samples). Must be numeric, character, or logical.")
    }
  } else {
    i <- seq_len(nsample(x))  # Default: all samples
  }
  
  # Validate and process 'j' indices (taxa)
  if (!missing(j)) {
    if (is.numeric(j)) {
      if (any(j < 1 | j > ntaxa(x))) {
        cli::cli_abort("{.var j} (taxa) indices out of range (1..{ntaxa(x)}).")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.character(j)) {
      j <- match(j, taxa_id(x))
      if (any(is.na(j))) {
        cli::cli_abort("Some taxa names in {.var j} not found in {.fn taxa_id}().")
      }
      check_unique(j, "j (taxa)")  # Check uniqueness
    } else if (is.logical(j)) {
      if (length(j) != ntaxa(x)) {
        cli::cli_abort("Logical {.var j} must have length equal to the number of taxa ({ntaxa(x)}).")
      }
      j <- which(j)  # Convert logical to numeric
    } else {
      cli::cli_abort("Invalid {.var j} (taxa). Must be numeric, character, or logical.")
    }
  } else {
    j <- seq_len(ntaxa(x))  # Default: all taxa
  }
  
  # Determine if full sample set is used to retain network and community
  full_i <- all(seq_len(nsample(x)) %in% i) || all(sample_id(x) %in% i)
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # Handle empty indices cases
  if ( length(i) == 0 && length(j) == 0 ){
    
    # return(new("mgnet")) 
    new_mgnet <- new("mgnet")
    
  } else if ( length(i) == 0 && length(j) > 0 ){
    
    # Case where only taxa are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(nrow=0, ncol=0)
    new_mgnet@rela <- matrix(nrow=0, ncol=0)
    new_mgnet@norm <- matrix(nrow=0, ncol=0)
    new_mgnet@meta <- data.frame()
    new_mgnet@taxa <- x@taxa[j,,drop = F]
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    # validObject(new_mgnet)
    # return(new_mgnet)
    
  } else if ( length(i) > 0 && length(j) == 0 ){
    
    # Case where only samples are selected
    new_mgnet <- x
    new_mgnet@abun <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@rela <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@norm <- matrix(numeric(0), nrow=0, ncol=0)
    new_mgnet@meta <- x@meta[i,,drop = F]
    new_mgnet@taxa <- data.frame()
    new_mgnet@netw <- igraph::make_empty_graph(0)
    new_mgnet@comm = igraph::cluster_fast_greedy(igraph::make_empty_graph(0,directed=F))
    # validObject(new_mgnet)
    # return(new_mgnet)
    
  } else {
    
    # Standard case where both i and j are present
    # Initialize new data subsets
    if(length(x@abun) != 0) abun.new <- x@abun[i,j,drop=F] else abun.new <- x@abun
    if(length(x@rela) != 0) rela.new <- x@rela[i,j,drop=F] else rela.new <- x@rela
    if(length(x@norm) != 0) norm.new <- x@norm[i,j,drop=F] else norm.new <- x@norm
    if(length(x@meta) != 0) meta.new <- x@meta[i, ,drop=F] else meta.new <- x@meta
    if(length(x@taxa) != 0) taxa.new <- x@taxa[j, ,drop=F] else taxa.new <- x@taxa
    
    if(length(x@netw) == 0){
      
      # When network is missing only abundances and metadata are returned
      new_mgnet <- mgnet(abun = abun.new,
                       rela = rela.new,
                       norm = norm.new,
                       meta = meta.new,
                       taxa = taxa.new)
      
    } else if (length(x@netw)!=0 & full_i) {
      
      netw <- x@netw
      netw.new <- igraph::induced_subgraph(netw, vids = V(netw)[j])
      netw.new <- igraph::permute(netw.new, match(igraph::V(netw.new)$name, igraph::V(netw)$name[j]))
      
      if(length(x@comm) != 0){
        v_old <- igraph::V(netw)$name
        v_new <- igraph::V(netw.new)$name
        
        m_new <- igraph::membership(x@comm)[match(v_new, v_old)]
        names(m_new) <- v_new
        
        algo <- attr(x@comm, "algorithm"); if (is.null(algo)) algo <- "subset"
        comm.new <- igraph::make_clusters(netw.new, membership = m_new,
                                          algorithm = algo, modularity = NA_real_)
      } else {
        comm.new <- x@comm
      }
      
      new_mgnet <- mgnet(abun = abun.new,
                       rela = rela.new,
                       norm = norm.new,
                       meta = meta.new,
                       taxa = taxa.new,
                       netw = netw.new,
                       comm = comm.new)
      
      if(are_selected_links(x)){
        old_selection <- get_selected_links(x)
        new_selection <- old_selection[old_selection %in% igraph::E(netw.new)$link_id]
        attr(new_mgnet, "selected_links") <- new_selection
      }
      
      if(is_link_grouped(x)){
        new_grouping <- get_grouped_link(x)[j]
        attr(new_mgnet, "link_groups") <- new_grouping
      }
      
      # return(new_mgnet)
      
    } else if (length(x@netw)!=0 & !full_i) {
      # SUBCASE NETWORK PRESENT AND I PRESENT
      
      cli::cli_warn("Subsetting by samples removes {.field netw} and {.field comm} slots.")
      new_mgnet <- mgnet(abun = abun.new,
                       rela = rela.new,
                       norm = norm.new,
                       meta = meta.new,
                       taxa = taxa.new)
      
    } else {
      
      cli::cli_abort("Why are you here? Only for suffering? (cit. Kaz).")
      
    }
    
  }

  if(is_mgnet_grouped(x)) new_mgnet <- group_mgnet(new_mgnet, !!!rlang::syms(get_group_mgnet(x)))
  if(is_link_grouped(x)) new_mgnet <- group_link(new_mgnet, !!!rlang::syms(get_grouped_link(x)))
  return(new_mgnet)
})