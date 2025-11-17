setMethod("show", "mgnet", function(object) {
  cat("==== mgnet Object Summary ====\n")
  
  # ------------ helpers ------------
  safe_int <- function(expr) { 
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) NA_integer_ else as.integer(out)
  }
  has_block <- function(x) !is.null(x) && length(x) > 0L &&
    is.matrix(x) && nrow(x) > 0L && ncol(x) > 0L
  
  # -------- available blocks -------
  A <- abun(object); R <- rela(object); N <- norm(object)
  present_blocks <- c(if (has_block(A)) "abun",
                      if (has_block(R)) "rela",
                      if (has_block(N)) "norm")
  
  # ---------- dimensions -----------
  ns <- safe_int(nsample(object))
  nt <- safe_int(ntaxa(object))
  cat("General Info:\n")
  cat(sprintf("  Samples: %s\n", ifelse(is.na(ns), "NA", ns)))
  cat(sprintf("  Taxa:    %s\n", ifelse(is.na(nt), "NA", nt)))
  
  # -- blocks statistics (priority abun>rela>norm) -------
  if (length(present_blocks)) {
    X <- if (has_block(A)) A else if (has_block(R)) R else N
    total_el <- as.double(nrow(X)) * as.double(ncol(X))
    if (total_el > 0) {
      nnz <- sum(X != 0, na.rm = TRUE)
      zero_pct <- 100 * (1 - (as.double(nnz) / total_el))
      na_cnt <- sum(is.na(X))
      cat(sprintf("  Active block: %s  |  Zeros: ~%.2f%%  |  NAs: %s\n",
                  paste(present_blocks, collapse = ","),
                  zero_pct, format(na_cnt, big.mark = ",")))
    } else {
      cat(sprintf("  Active block: %s  |  Zeros: NA%%  |  NAs: NA\n",
                  paste(present_blocks, collapse = ",")))
    }
  } else {
    cat("  No abundance-like data available (abun/rela/norm all empty).\n")
  }
  
  # ------------ meta / taxa ------------
  M <- meta(object)
  if (!is.null(M) && is.data.frame(M) && ncol(M) > 0L) {
    cols <- utils::head(colnames(M), 6)
    cat(sprintf("  Sample meta cols: %s%s\n",
                toString(cols), if (ncol(M) > 6) ", ..." else ""))
  } else cat("  Sample meta: none.\n")
  
  Tt <- taxa(object)
  if (!is.null(Tt) && is.data.frame(Tt) && ncol(Tt) > 0L) {
    cols <- utils::head(colnames(Tt), 6)
    cat(sprintf("  Taxa meta cols:   %s%s\n",
                toString(cols), if (ncol(Tt) > 6) ", ..." else ""))
  } else cat("  Taxa meta: none.\n")
  
  # ------------ network & selected links ------------
  rawG <- object@netw
  raw_nodes <- if (inherits(rawG, "igraph")) igraph::vcount(rawG) else 0L
  raw_edges <- if (inherits(rawG, "igraph")) igraph::ecount(rawG) else 0L
  
  sel_active <- isTRUE(are_selected_links(object))
  sel_ids    <- if (sel_active) get_selected_links(object) else integer(0)
  selE       <- length(sel_ids)
  sel_pct    <- if (sel_active && raw_edges > 0L) 100 * selE / raw_edges else NA_real_
  
  G <- netw(object)
  if (inherits(G, "igraph") && igraph::vcount(G) > 0L) {
    dens_shown <- 100 * igraph::edge_density(G)
    cat(sprintf("  Network (shown): %d nodes, %d edges, density ~%.2f%%\n",
                igraph::vcount(G), igraph::ecount(G), dens_shown))
    if (sel_active && raw_edges > 0L) {
      cat(sprintf("    Links selected: YES (%d of %d, ~%.2f%%). Showing selected subgraph.\n",
                  selE, raw_edges, sel_pct))
    } else if (raw_edges > 0L) {
      cat("    Links selected: NO (showing full graph).\n")
    }
  } else {
    if (raw_nodes > 0L) {
      cat(sprintf("  Network: empty after filtering. Raw graph had %d nodes, %d edges.\n",
                  raw_nodes, raw_edges))
      if (sel_active && raw_edges > 0L) {
        cat(sprintf("    Links selected: YES (%d of %d, ~%.2f%%). Selection produced empty view.\n",
                    selE, raw_edges, sel_pct))
      }
    } else {
      cat("  Network: none.\n")
    }
  }
  
  # ------------ communities ------------
  Cc <- comm(object)
  if (!is.null(Cc)) {
    sz <- igraph::sizes(Cc)
    if (length(sz)) {
      sz <- as.integer(sz)
      n_iso <- sum(sz == 1L)
      n_comm_gt1 <- sum(sz > 1L)
      shown_nodes <- if (inherits(G, "igraph")) igraph::vcount(G) else 0L
      iso_pct <- if (shown_nodes > 0L) 100 * (n_iso / shown_nodes) else NA_real_
      if (n_comm_gt1 > 0L) {
        top_sz <- utils::head(sort(sz[sz > 1L], decreasing = TRUE), 6)
        cat(sprintf("  Communities (>1): %d  |  Isolated: %d (~%s)\n",
                    n_comm_gt1, n_iso,
                    ifelse(is.na(iso_pct), "NA%", sprintf("%.2f%%", iso_pct))))
        if (length(top_sz)) {
          cat(sprintf("  Top community sizes: %s%s\n",
                      toString(top_sz),
                      if (sum(sz > 1L) > length(top_sz)) ", ..." else ""))
        }
      } else {
        cat("  Communities: all nodes isolated.\n")
      }
    } else cat("  Communities: empty.\n")
  } else cat("  Communities: none.\n")
  
  cat("==== End of mgnet Object Summary ====\n")
})

setMethod("show", "mgnets", function(object) {
  cat("==== mgnets Object Summary ====\n")
  k  <- length(object)
  nm <- names(object); has_names <- !is.null(nm) && length(nm) == k
  cat(sprintf("Elements: %d\n", k))
  if (k == 0L) {
    cat("Empty mgnets: no mgnet objects.\n")
    cat("==== End of mgnets Object Summary ====\n")
    return(invisible(NULL))
  }
  
  safe_int <- function(expr) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) NA_integer_ else as.integer(out)
  }
  has_block <- function(x) !is.null(x) && length(x) > 0L &&
    is.matrix(x) && nrow(x) > 0L && ncol(x) > 0L
  
  cat("Preview (up to 5 elements):\n")
  for (i in seq_len(min(k, 5L))) {
    tag <- if (has_names && nzchar(nm[i])) sprintf("'%s'", nm[i]) else sprintf("#%d", i)
    om  <- object[[i]]
    
    if (!methods::is(om, "mgnet")) {
      cat(sprintf("  - %s: [not an 'mgnet' instance]\n", tag)); next
    }
    
    ns <- safe_int(nsample(om))
    nt <- safe_int(ntaxa(om))
    
    A <- abun(om); R <- rela(om); N <- norm(om)
    present_blocks <- c(if (has_block(A)) "abun",
                        if (has_block(R)) "rela",
                        if (has_block(N)) "norm")
    X <- if (has_block(A)) A else if (has_block(R)) R else if (has_block(N)) N else NULL
    
    zero_pct <- NA_real_; na_cnt <- NA_integer_
    if (!is.null(X)) {
      total_el <- as.double(nrow(X)) * as.double(ncol(X))
      if (total_el > 0) {
        nnz <- sum(X != 0, na.rm = TRUE)
        zero_pct <- 100 * (1 - (as.double(nnz) / total_el))
        na_cnt <- sum(is.na(X))
      }
    }
    
    rawG <- om@netw
    raw_edges <- if (inherits(rawG, "igraph")) igraph::ecount(rawG) else 0L
    sel_active <- isTRUE(are_selected_links(om))
    sel_ids    <- if (sel_active) get_selected_links(om) else integer(0)
    selE       <- length(sel_ids)
    sel_pct    <- if (sel_active && raw_edges > 0L) 100 * selE / raw_edges else NA_real_
    
    G <- netw(om)
    shown_nodes <- if (inherits(G, "igraph")) igraph::vcount(G) else 0L
    shown_edges <- if (inherits(G, "igraph")) igraph::ecount(G) else 0L
    dens_shown  <- if (inherits(G, "igraph") && shown_nodes > 0L) 100 * igraph::edge_density(G) else NA_real_
    
    Cc <- comm(om)
    n_comm_gt1 <- NA_integer_; iso_pct <- NA_real_
    if (!is.null(Cc)) {
      sz <- igraph::sizes(Cc)
      if (length(sz)) {
        sz <- as.integer(sz)
        n_comm_gt1 <- sum(sz > 1L)
        if (shown_nodes > 0L) iso_pct <- 100 * (sum(sz == 1L) / shown_nodes)
      }
    }
    
    cat(sprintf("  - %s\n", tag))
    cat(sprintf("      dims: samples=%s, taxa=%s\n",
                ifelse(is.na(ns), "NA", ns), ifelse(is.na(nt), "NA", nt)))
    cat(sprintf("      data: active=%s%s\n",
                if (length(present_blocks)) paste(present_blocks, collapse = ",") else "none",
                if (!is.na(zero_pct) || !is.na(na_cnt))
                  sprintf(" | zeros~%s%s",
                          ifelse(is.na(zero_pct), "NA%",
                                 sprintf("%.2f%%", zero_pct)),
                          if (!is.na(na_cnt)) sprintf(", NAs=%s", format(na_cnt, big.mark=",")) else "")
                else ""))
    link_tag <- if (sel_active && raw_edges > 0L) {
      sprintf("selected %d/%d (~%.2f%%)", selE, raw_edges, sel_pct)
    } else if (raw_edges > 0L) "all" else "none"
    cat(sprintf("      netw: nodes=%d, edges=%d, dens~%s | links=%s\n",
                shown_nodes, shown_edges,
                ifelse(is.na(dens_shown), "NA%", sprintf("%.2f%%", dens_shown)),
                link_tag))
    cat(sprintf("      comm: groups>1=%s, iso~%s\n",
                ifelse(is.na(n_comm_gt1), "NA", n_comm_gt1),
                ifelse(is.na(iso_pct), "NA%", sprintf("%.2f%%", iso_pct))))
  }
  
  if (k > 5L) cat(sprintf("  ... and %d more.\n", k - 5L))
  cat("==== End of mgnets Object Summary ====\n")
})


