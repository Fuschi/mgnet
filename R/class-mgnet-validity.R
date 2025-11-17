################################################################################
##                             SHARED SLOT CHECKS                             ##
################################################################################
#' @include class-mgnet.R
NULL

#------------------------------------------------------------------------------#
#' Assert: unique row and column names
#' @keywords internal
.assert_named_and_uniqueness <- function(x, field) {
  rn <- rownames(x); cn <- colnames(x)
  
  if (is.null(rn)) cli::cli_abort("The {.field {field}} must have row names.", class="mgnet_validators_error")
  if (anyNA(rn) || any(!nzchar(rn))) cli::cli_abort("Row names in {.field {field}} must be non-empty and non-NA.", class="mgnet_validators_error")
  if (anyDuplicated(rn) > 0) cli::cli_abort("Duplicate row names in {.field {field}}.", class="mgnet_validators_error")
  
  if (is.null(cn)) cli::cli_abort("The {.field {field}} must have column names.", class="mgnet_validators_error")
  if (anyNA(cn) || any(!nzchar(cn))) cli::cli_abort("Column names in {.field {field}} must be non-empty and non-NA.", class="mgnet_validators_error")
  if (anyDuplicated(cn) > 0) cli::cli_abort("Duplicate column names in {.field {field}}.", class="mgnet_validators_error")
}

#------------------------------------------------------------------------------#
#' Assert: no reserved keywords in column names
#' @keywords internal
.assert_no_reserved_names <- function(x, field) {
  conflicted <- intersect(colnames(x), .MGNET_RESERVED_COLNAMES)
  if (length(conflicted) > 0) {
    cli::cli_abort(
      c(
        "Column name{?s} in {.field {field}} must not match reserved keyword{?s}.",
        "x" = "Found: {.val {conflicted}}"
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: only values >=0
#' @keywords internal
.assert_positive_values <- function(x, field) {
  if (any(x < 0, na.rm = TRUE)) {
    cli::cli_abort("Negative values found in {.field {field}}; all values must be >= 0.",
                   class = "mgnet_validators_error")
  }
}

################################################################################
##                             SHARED RECIPROCAL CHECKS                       ##
################################################################################

#------------------------------------------------------------------------------#
#' Assert: reciprocal names of rows and columns must be identical
#' @keywords internal
.assert_reciprocal_dimensions <- function(x1, field1, x2, field2) {
  if (nrow(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Row count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Rows represent samples. Ensure both have the same number of samples.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (ncol(x1) != ncol(x2)) {
    cli::cli_abort(
      c(
        "Column count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Columns represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field {field2}} has {.val {ncol(x2)}} columns."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(rownames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Row names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Row names correspond to sample IDs and must match exactly in content and order."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(colnames(x1), colnames(x2))) {
    cli::cli_abort(
      c(
        "Column names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Column names correspond to taxa IDs and must match exactly in content and order."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of rows
#' @keywords internal
.assert_identical_samples <- function(x1, field1, x2, field2) {
  if (nrow(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Row count mismatch between {.field {field1}} and {.field {field2}}.",
        "i" = "Rows represent samples. Ensure both have the same number of samples.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(rownames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Row names of {.field {field1}} and {.field {field2}} must be identical.",
        "i" = "Row names correspond to sample IDs and must match exactly in content and order."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of cols of x1 and rows of x2
#' @keywords internal
.assert_identical_taxa <- function(x1, field1, x2, field2) {
  if (ncol(x1) != nrow(x2)) {
    cli::cli_abort(
      c(
        "Columns count in {.field {field1}} and rows in {.field {field2}} mismatch.",
        "i" = "Columns in {.field {field1}}, rows in {.field {field2}} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field {field2}} has {.val {nrow(x2)}} rows."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(colnames(x1), rownames(x2))) {
    cli::cli_abort(
      c(
        "Column names in {.field {field1}} and row names in {.field {field2}} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of cols of x1 and nodes in netw
#' @keywords internal
.assert_identical_abun_nodes <- function(x1, field1, netw) {
  if (ncol(x1) != igraph::vcount(netw)) {
    cli::cli_abort(
      c(
        "Columns count in {.field {field1}} and nodes in {.field netw} mismatch.",
        "i" = "Columns in {.field {field1}}, nodes in {.field netw} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {ncol(x1)}} columns.",
        "x" = "{.field netw} has {.val {igraph::vcount(netw)}} nodes."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(colnames(x1), igraph::V(netw)$name)) {
    cli::cli_abort(
      c(
        "Column names in {.field {field1}} and nodes names in {.field netw} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: reciprocal number and names of rows of x1 and nodes in netw
#' @keywords internal
.assert_identical_taxa_nodes <- function(x1, field1, netw) {
  if (nrow(x1) != igraph::vcount(netw)) {
    cli::cli_abort(
      c(
        "Rows count in {.field {field1}} and nodes in {.field netw} mismatch.",
        "i" = "Rows in {.field {field1}}, nodes in {.field netw} represent taxa. Ensure both have the same number of taxa.",
        "x" = "{.field {field1}} has {.val {nrow(x1)}} rows.",
        "x" = "{.field netw} has {.val {igraph::vcount(netw)}} nodes."
      ),
      class = "mgnet_validators_error"
    )
  }
  if (!identical(rownames(x1), igraph::V(netw)$name)) {
    cli::cli_abort(
      c(
        "Row names in {.field {field1}} and node names in {.field netw} must match.",
        "i" = "These names represent taxa and must be identical and in the same order."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: identical zero position in the two matrices
#' @keywords internal
.assert_abun_zero_structure <- function(x1, field1, x2, field2) {
  if (!identical(x1 == 0, x2 == 0)) {
    cli::cli_abort(
      c(
        "Inconsistent zero patterns between {.field {field1}} and {.field {field2}}.",
        "i" = "Zeros must appear in the same positions in both matrices."
      ),
      class = "mgnet_validators_error"
    )
  }
}

#------------------------------------------------------------------------------#
#' Assert: consistency between network and community
#' @keywords internal
.assert_netw_comm_consistency <- function(netw, comm) {
  if (length(netw) == 0 && length(comm) > 0) {
    cli::cli_abort("The {.field comm} slot cannot be present without a {.field netw} slot.",
                   class = "mgnet_validators_error")
  }
  if (length(netw) > 0 && length(comm) > 0) {
    if (length(comm$membership) != igraph::vcount(netw)) {
      cli::cli_abort("The number of membership elements in {.field comm} must match the number of vertices in {.field netw}.",
                     class = "mgnet_validators_error")
    }
  }
}

#------------------------------------------------------------------------------#
#' Assert: disjoint feature names between meta and taxa
#' @keywords internal
.assert_disjoint_feature_names <- function(meta, taxa) {
  cnm <- colnames(meta)
  cnt <- colnames(taxa)
  conflicted <- intersect(cnm, cnt)
  if (length(conflicted) > 0) {
    cli::cli_abort(
      c(
        "Column names of {.field meta} and {.field taxa} must be disjoint.",
        "x" = "Found in both: {.val {conflicted}}"
      ),
      class = "mgnet_validators_error"
    )
  }
}

################################################################################
##                                SLOT CHECKS                                 ##
################################################################################

#------------------------------------------------------------------------------#
#' Validate abun slot
#' @keywords internal
.validate_abun <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field abun} slot must be a numeric matrix (samples x taxa).",
                   class = "mgnet_validators_error")
  }
  .assert_named_and_uniqueness(x, "abun")
  .assert_no_reserved_names(x, "abun")
  .assert_positive_values(x, "abun")
}

#------------------------------------------------------------------------------#
#' Validate rela slot
#' @keywords internal
.validate_rela <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field rela} slot must be a numeric matrix (samples x taxa).",
                   class = "mgnet_validators_error")
  }
  .assert_named_and_uniqueness(x, "rela")
  .assert_no_reserved_names(x, "rela")
  .assert_positive_values(x, "rela")
}

#------------------------------------------------------------------------------#
#' Validate norm slot
#' @keywords internal
.validate_norm <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("The {.field norm} slot must be a numeric matrix (samples x taxa).",
                   class = "mgnet_validators_error")
  }
  .assert_named_and_uniqueness(x, "norm")
  .assert_no_reserved_names(x, "norm")
}

#------------------------------------------------------------------------------#
#' Validate meta slot
#' @keywords internal
.validate_meta <- function(x) {
  if (!is.data.frame(x)) {
    cli::cli_abort("The {.field meta} slot must be a {.cls data.frame} with sample metadata (samples x features).",
                   class = "mgnet_validators_error")
  }
  .assert_named_and_uniqueness(x, "meta")
  .assert_no_reserved_names(x, "meta")
}

#------------------------------------------------------------------------------#
#' Validate taxa slot
#' @keywords internal
.validate_taxa <- function(x) {
  if (!is.data.frame(x)) {
    cli::cli_abort("The {.field taxa} slot must be a {.cls data.frame} with taxa metadata (taxa x features).",
                   class = "mgnet_validators_error")
  }
  .assert_named_and_uniqueness(x, "taxa")
  .assert_no_reserved_names(x, "taxa")
}

#------------------------------------------------------------------------------#
#' Validate netw slot
#' @keywords internal
.validate_netw <- function(x) {
  if (!igraph::is_igraph(x)) {
    cli::cli_abort("The {.field netw} slot must be an {.cls igraph}.", class="mgnet_validators_error")
  }
  if (!igraph::is_named(x)) {
    cli::cli_abort("The {.field netw} graph must have vertex names.", class="mgnet_validators_error")
  }
  vnames <- igraph::V(x)$name
  if (anyNA(vnames) || any(!nzchar(vnames)) || anyDuplicated(vnames) > 0) {
    cli::cli_abort("Vertex names in {.field netw} must be unique, non-empty and non-NA.", class="mgnet_validators_error")
  }
  
  # -- link_id checks ----------------------------------------------------------
  m <- igraph::ecount(x)
  if (m == 0L) return(invisible())
  
  lid <- igraph::edge_attr(x, "link_id")
  if (is.null(lid)) {
    cli::cli_abort(
      "Edge attribute {.val link_id} is missing in {.field netw}.",
      class = "mgnet_validators_error"
    )
  }
  
  if (length(lid) != m) {
    cli::cli_abort(
      "Edge attribute {.val link_id} must have length {m}, found {length(lid)}.",
      class = "mgnet_validators_error"
    )
  }
  
  if (!(is.numeric(lid) || is.character(lid))) {
    cli::cli_abort(
      "{.val link_id} must be either numeric/integer or character.",
      class = "mgnet_validators_error"
    )
  }
  
  if (anyNA(lid)) {
    cli::cli_abort(
      "{.val link_id} cannot contain NA values.",
      class = "mgnet_validators_error"
    )
  }
  
  if (is.character(lid) && any(!nzchar(lid))) {
    cli::cli_abort(
      "{.val link_id} cannot contain empty strings.",
      class = "mgnet_validators_error"
    )
  }
  
  if (anyDuplicated(lid) > 0) {
    cli::cli_abort(
      "{.val link_id} must be unique for each edge.",
      class = "mgnet_validators_error"
    )
  }

}

#------------------------------------------------------------------------------#
#' Validate comm slot
#' @keywords internal
.validate_comm <- function(x) {
  if (!inherits(x, "communities")) {
    cli::cli_abort("The {.field comm} slot must be a {.cls communities} object.",
                   class = "mgnet_validators_error")
  }
}

################################################################################
##                             CHECK RECIPROCAL                               ##
################################################################################
#' Validate comm slot
#' @keywords internal
.validate_reciprocal <- function(obj) {
  
  if (length(obj@abun) > 0){
    if (length(obj@rela) > 0) {.assert_reciprocal_dimensions(obj@abun, "abun", obj@rela, "rela")
                               .assert_abun_zero_structure(obj@abun, "abun", obj@rela, "rela")}
    if(length(obj@norm) > 0) .assert_reciprocal_dimensions(obj@abun, "abun", obj@norm, "norm")
    if(length(obj@meta) > 0) .assert_identical_samples(obj@abun, "abun", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@abun, "abun", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@abun, "abun", obj@netw)
  } 
  
  if (length(obj@rela) > 0){
    if(length(obj@norm) > 0) .assert_reciprocal_dimensions(obj@rela, "rela", obj@norm, "norm")
    if(length(obj@meta) > 0) .assert_identical_samples(obj@rela, "rela", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@rela, "rela", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@rela, "rela", obj@netw)
  } 
  
  if (length(obj@norm) > 0){
    if(length(obj@meta) > 0) .assert_identical_samples(obj@norm, "norm", obj@meta, "meta")
    if(length(obj@taxa) > 0) .assert_identical_taxa(obj@norm, "norm", obj@taxa, "taxa")
    if(length(obj@netw) > 0) .assert_identical_abun_nodes(obj@norm, "norm", obj@netw)
  } 
  
  if (length(obj@taxa) > 0){
    if(length(obj@netw) > 0) .assert_identical_taxa_nodes(obj@taxa, "taxa", obj@netw)
  } 
  
  if (length(obj@meta) > 0 && length(obj@taxa) > 0) {
    .assert_disjoint_feature_names(obj@meta, obj@taxa)
  }
  
  .assert_netw_comm_consistency(obj@netw, obj@comm)
  
}

################################################################################
##                           GROUPING/SELECTION                               ##
################################################################################

#------------------------------------------------------------------------------#
#' Validate mgnet_groups attribute (grouping on meta/taxa)
#' @keywords internal
.validate_mgnet_groups <- function(object) {
  groups <- attr(object, "mgnet_groups", exact = TRUE)
  
  # No groups: nothing to control
  if (is.null(groups)) return(invisible())
  
  # Type
  if (!is.character(groups)) {
    cli::cli_abort(
      c(
        "x" = 'The "mgnet_groups" attribute must be a character vector.',
        "i" = "Use {.fn group_mgnet}() to set grouping variables."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  # NA / empty strings
  if (anyNA(groups) || any(!nzchar(groups))) {
    cli::cli_abort(
      'Grouping variables in "mgnet_groups" must be non-NA, non-empty strings.',
      class = "mgnet_validators_error"
    )
  }
  
  # Duplicates
  if (anyDuplicated(groups) > 0L) {
    cli::cli_abort(
      'Grouping variables in "mgnet_groups" must be unique.',
      class = "mgnet_validators_error"
    )
  }
  
  # Coherence with current meta/taxa variables
  available <- meta_taxa_vars(object)
  missing   <- setdiff(groups, available)
  
  if (length(missing)) {
    cli::cli_abort(
      c(
        "x" = 'Unknown grouping variable{?s} stored in attribute {.field "mgnet_groups"}: {.val {missing}}.',
        "v" = "Available variables are: {.val {available}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  invisible()
}

#------------------------------------------------------------------------------#
#' Validate selected_links attribute (edge selection)
#' @keywords internal
.validate_selected_links <- function(object) {
  sel <- attr(object, "selected_links", exact = TRUE)
  
  # No selection: nothing to check
  if (is.null(sel)) return(invisible())
  
  # If there is a selection but no network -> error
  if (miss_netw(object)) {
    cli::cli_abort(
      c(
        "x" = 'Attribute {.field "selected_links"} is set, but no network is present in {.cls mgnet}.',
        "i" = "Either add a network or call {.fn deselect_link}() before validating."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  # Must be an atomic vector (typically numeric or character)
  if (!is.atomic(sel)) {
    cli::cli_abort(
      c(
        "x" = 'Attribute {.field "selected_links"} must be an atomic vector of link IDs.',
        "v" = "Got: {.cls {class(sel)}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  if (anyNA(sel)) {
    cli::cli_abort(
      'Attribute {.field "selected_links"} must not contain NA values.',
      class = "mgnet_validators_error"
    )
  }
  
  # Retrieve current edges and their link_id
  edges <- igraph::as_data_frame(object@netw, what = "edges")
  
  if (!"link_id" %in% names(edges)) {
    cli::cli_abort(
      c(
        "x" = 'Network edges must contain a {.field "link_id"} column to support edge selection.',
        "i" = "Ensure that {.fn .link_prepare}() or the network constructor sets a link_id edge attribute."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  valid_ids <- edges[["link_id"]]
  
  # Optional: coerce to same type if needed (e.g. both to character)
  sel_coerced    <- as.character(sel)
  valid_ids_comp <- as.character(valid_ids)
  missing_ids    <- setdiff(sel_coerced, valid_ids_comp)
  
  missing_ids <- setdiff(sel, valid_ids)
  
  if (length(missing_ids)) {
    cli::cli_abort(
      c(
        "x" = 'Some values in {.field "selected_links"} do not match any current edge {.field "link_id"}.',
        "v" = "Unknown IDs: {.val {missing_ids}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  invisible()
}


#------------------------------------------------------------------------------#
#' Validate link_groups attribute (edge-level grouping)
#' @keywords internal
.validate_link_groups <- function(object) {
  groups <- attr(object, "link_groups", exact = TRUE)
  
  # No groups: nothing to check
  if (is.null(groups)) return(invisible())
  
  # If I have a grouping without a network -> error
  if (miss_netw(object)) {
    cli::cli_abort(
      c(
        "x" = 'Attribute {.field "link_groups"} is set, but no network is present in {.cls mgnet}.',
        "i" = "Either add a network or call {.fn ungroup_link}() before validating."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  if (!is.numeric(groups)) {
    cli::cli_abort(
      c(
        "x" = 'Attribute {.field "link_groups"} must be a numeric (integer) vector of group IDs.',
        "v" = "Got: {.cls {class(groups)}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  if (anyNA(groups) || any(!is.finite(groups))) {
    cli::cli_abort(
      'Attribute {.field "link_groups"} must not contain NA or non-finite values.',
      class = "mgnet_validators_error"
    )
  }
  
  if (any(groups < 1)) {
    cli::cli_abort(
      'Attribute {.field "link_groups"} must contain strictly positive group IDs.',
      class = "mgnet_validators_error"
    )
  }
  
  # The groups must be of the same length of the number of the edges.
  g <- netw(object, selected = FALSE)
  m <- igraph::ecount(g)
  
  if (length(groups) != m) {
    cli::cli_abort(
      c(
        "x" = 'Length of {.field "link_groups"} must match the number of edges in the network.',
        "v" = "{.field \"link_groups\"} length: {.val {length(groups)}}, edges in network: {.val {m}}."
      ),
      class = "mgnet_validators_error"
    )
  }
  
  invisible()
}


################################################################################
##                                 VALIDATOR                                  ##
################################################################################

#------------------------------------------------------------------------------#
#' Validity method for mgnet class
#' @name mgnet-validator
#' @keywords internal
setValidity("mgnet", function(object) {
  tryCatch({
    if (length(object@abun) > 0) .validate_abun(object@abun)
    if (length(object@rela) > 0) .validate_rela(object@rela)
    if (length(object@norm) > 0) .validate_norm(object@norm)
    if (length(object@meta) > 0) .validate_meta(object@meta)
    if (length(object@taxa) > 0) .validate_taxa(object@taxa)
    if (length(object@netw) > 0) .validate_netw(object@netw)
    if (length(object@comm) > 0) .validate_comm(object@comm)
    .validate_reciprocal(object)
    .validate_mgnet_groups(object)
    .validate_selected_links(object)
    .validate_link_groups(object)
    TRUE
  }, mgnet_validators_error = function(e) {
    conditionMessage(e)
  })
})
#------------------------------------------------------------------------------#
# End of validity definition


