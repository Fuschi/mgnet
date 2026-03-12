#'@include class-mgnet.R class-mgnets.R
NULL

# Internal helpers shared across setters
#------------------------------------------------------------------------------#
.normalize_input_ids <- function(value, id_col, allow_matrix = TRUE) {
  is_tabular <- inherits(value, "data.frame") || inherits(value, "tbl_df")
  is_mat <- inherits(value, "matrix")
  
  if (!(is_tabular || (allow_matrix && is_mat))) {
    cli::cli_abort("{.arg value} must be a matrix, tibble, or data frame.")
  }
  
  has_id_col <- id_col %in% colnames(value)
  has_rownames <- if (is_tabular) tibble::has_rownames(value) else !is.null(rownames(value))
  
  if (has_id_col && has_rownames) {
    cli::cli_abort(
      "{.arg value} cannot have both a {.val {id_col}} column and row names."
    )
  }
  
  if (!has_id_col && !has_rownames) {
    cli::cli_abort(
      "{.arg value} must have either a {.val {id_col}} column or row names."
    )
  }
  
  if (has_id_col) {
    ids <- value[[id_col]]
    
    if (anyDuplicated(ids)) {
      cli::cli_abort(
        "Column {.var {id_col}} contains duplicated values in {.arg value}."
      )
    }
    
    if (is_tabular) {
      value <- tibble::column_to_rownames(as.data.frame(value), id_col)
    } else {
      rownames(value) <- ids
      value <- value[, colnames(value) != id_col, drop = FALSE]
    }
  }
  
  value
}

.align_to_ids <- function(value, expected_ids, id_col) {
  current_ids <- rownames(value)
  
  if (!all(current_ids %in% expected_ids)) {
    cli::cli_abort(
      "Provided {.var {id_col}} values do not match those in the object."
    )
  }
  
  if (!identical(current_ids, expected_ids)) {
    value <- value[expected_ids, , drop = FALSE]
  }
  
  value
}

.assign_slot_checked <- function(object, slot_name, value, coercer) {
  methods::slot(object, slot_name) <- coercer(value)
  methods::validObject(object)
  object
}

.set_sample_matrix_slot <- function(object, value, slot_name) {
  value <- .normalize_input_ids(value, id_col = "sample_id", allow_matrix = TRUE)
  
  expected_ids <- sample_id(object)
  if (length(expected_ids) > 0L) {
    value <- .align_to_ids(value, expected_ids, "sample_id")
  }
  
  .assign_slot_checked(object, slot_name, value, as.matrix)
}

.set_annotation_slot <- function(object, value, slot_name, id_col, expected_ids = NULL) {
  value <- .normalize_input_ids(value, id_col = id_col, allow_matrix = FALSE)
  
  if (!is.null(expected_ids) && length(expected_ids) > 0L) {
    value <- .align_to_ids(value, expected_ids, id_col)
  }
  
  .assign_slot_checked(object, slot_name, value, as.data.frame)
}

.apply_mgnets_setter_list <- function(object, value, setter_name) {
  is_list_mgnets_assign(object, value)
  
  setter <- get(setter_name, mode = "function")
  for (nm in names(object)) {
    object@mgnets[[nm]] <- setter(object@mgnets[[nm]], value[[nm]])
  }
  
  methods::validObject(object)
  object
}

.split_mgnets_table_value <- function(object, value, level = c("sample", "taxa")) {
  level <- match.arg(level)
  
  if (is.list(value) && !is.data.frame(value)) {
    is_list_mgnets_assign(object, value)
    return(value)
  }
  
  if (!is.data.frame(value)) {
    cli::cli_abort(
      "{.arg value} must be either a named list of data frames or a single data frame."
    )
  }
  
  is_assign_mgnets_tbl(object, value, level)
  pieces <- split(value, value$mgnet)
  
  lapply(pieces, function(x) {
    x$mgnet <- NULL
    x
  })
}

.set_mgnets_annotation_slot <- function(object, value, setter_name, level = c("sample", "taxa")) {
  pieces <- .split_mgnets_table_value(object, value, level = level)
  setter <- get(setter_name, mode = "function")
  
  for (nm in names(object)) {
    object@mgnets[[nm]] <- setter(object@mgnets[[nm]], pieces[[nm]])
  }
  
  methods::validObject(object)
  object
}

# ABUNDANCE
#------------------------------------------------------------------------------#
#' Set Abundance Data
#'
#' This setter function allows you to update the abundance data for `mgnet` objects
#' and each `mgnet` object within a `mgnets`. The abundance data must be a numeric matrix
#' for `mgnet` objects. For `mgnets` objects, the abundance data should be a named list
#' of numeric matrices corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param value The new abundance data to be set.
#' @return The `mgnet` or `mgnets` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name abun<-
#' @aliases abun<-,mgnet-method abun<-,mgnets-method
setGeneric("abun<-", function(object, value) standardGeneric("abun<-"))

setMethod("abun<-", c("mgnet","ANY"), function(object, value){
  .set_sample_matrix_slot(object, value, "abun")
})

setMethod("abun<-", c("mgnets","ANY"), function(object, value){
  .apply_mgnets_setter_list(object, value, "abun<-")
})


# rela
#------------------------------------------------------------------------------#
#' Set rela Data
#'
#' This setter function allows you to update the rela data for `mgnet` objects
#' and each `mgnet` object within a `mgnets`. The rela data must be a numeric matrix
#' for `mgnet` objects. For `mgnets` objects, the rela data should be a named list
#' of numeric matrices corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param value The new rela data to be set.
#' @return The `mgnet` or `mgnets` object with updated abundance data.
#' @export
#' @importFrom methods validObject
#' @name rela<-
#' @aliases rela<-,mgnet-method rela<-,mgnets-method
setGeneric("rela<-", function(object, value) standardGeneric("rela<-"))

setMethod("rela<-", c("mgnet","ANY"), function(object, value){
  .set_sample_matrix_slot(object, value, "rela")
})

setMethod("rela<-", c("mgnets","ANY"), function(object, value){
  .apply_mgnets_setter_list(object, value, "rela<-")
})


# norm
#------------------------------------------------------------------------------#
#' Set Log-Transformed Abundance Data
#'
#' This setter function allows you to update the log-transformed abundance data for `mgnet` objects
#' and each `mgnet` object within a `mgnets`. The log-transformed abundance data must be a numeric matrix
#' for `mgnet` objects. For `mgnets` objects, the data should be a named list of numeric matrices
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param value The new log-transformed abundance data to be set.
#' @return The `mgnet` or `mgnets` object with updated log-transformed abundance data.
#' @export
#' @importFrom methods validObject
#' @name norm<-
#' @aliases norm<-,mgnet-method norm<-,mgnets-method
setGeneric("norm<-", function(object, value) standardGeneric("norm<-"))

setMethod("norm<-", c("mgnet", "ANY"), function(object, value) {
  .set_sample_matrix_slot(object, value, "norm")
})

setMethod("norm<-", c("mgnets", "ANY"), function(object, value) {
  .apply_mgnets_setter_list(object, value, "norm<-")
})



# META
#------------------------------------------------------------------------------#
#' Update Sample Metadata for mgnet and mgnets Objects
#'
#' This function sets the sample metadata for `mgnet` objects and also updates
#' each `mgnet` object within a `mgnets`. The metadata for `mgnet` objects must be
#' provided as a dataframe. For `mgnets` objects, the metadata should be either
#' a named list of dataframes corresponding to each `mgnet` object within the list,
#' or a single dataframe containing a 'mgnet' and 'sample_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `mgnet` or `mgnets` object to be updated.
#' @param value The new sample metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `mgnet` or `mgnets` object with new sample metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames has_rownames
#' @name meta<-
#' @aliases meta<-,mgnet-method meta<-,mgnets-method
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))

setMethod("meta<-", c("mgnet", "ANY"), function(object, value) {
  expected_ids <- sample_id(object)
  if (length(expected_ids) == 0L) expected_ids <- NULL
  
  .set_annotation_slot(
    object = object,
    value = value,
    slot_name = "meta",
    id_col = "sample_id",
    expected_ids = expected_ids
  )
})

setMethod("meta<-", c("mgnets","ANY"), function(object, value){
  .set_mgnets_annotation_slot(object, value, "meta<-", level = "sample")
})


# TAXA
#------------------------------------------------------------------------------#
#' Update Taxa Metadata for mgnet and mgnets Objects
#'
#' This function sets the taxa metadata for `mgnet` objects and also updates
#' each `mgnet` object within a `mgnets`. The metadata for `mgnet` objects must be
#' provided as a dataframe. For `mgnets` objects, the metadata should be either
#' a named list of dataframes corresponding to each `mgnet` object within the list,
#' or a single dataframe containing a 'mgnet' and 'taxa_id' columns to split the 
#' metadata accordingly.
#'
#' @param object An `mgnet` or `mgnets` object to be updated.
#' @param value The new taxa metadata to set, either a dataframe or a named list
#' of dataframes as per the object type.
#' @return The updated `mgnet` or `mgnets` object with new taxa metadata.
#' @export
#' @importFrom methods validObject
#' @importFrom tibble column_to_rownames
#' @name taxa<-
#' @aliases taxa<-,mgnet-method taxa<-,mgnets-method
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))

setMethod("taxa<-", c("mgnet", "ANY"), function(object, value) {
  expected_ids <- taxa_id(object)
  if (length(expected_ids) == 0L) expected_ids <- NULL
  
  value <- .normalize_input_ids(value, id_col = "taxa_id", allow_matrix = FALSE)
  
  if (!is.null(expected_ids)) {
    value <- .align_to_ids(value, expected_ids, "taxa_id")
  }
  
  if ("comm_id" %in% colnames(value)) {
    if (!all(value$comm_id == comm_id(object))) {
      cli::cli_abort("Community memberships in column {.var comm_id} differ from those stored in the object.")
    }
    value$comm_id <- NULL
  }
  
  .assign_slot_checked(object, "taxa", value, as.data.frame)
})

setMethod("taxa<-", c("mgnets","ANY"), function(object, value){
  .set_mgnets_annotation_slot(object, value, "taxa<-", level = "taxa")
})


# NETWORK
#------------------------------------------------------------------------------#
#' Set Network Data
#'
#' This setter function updates the network data for `mgnet` objects and
#' each `mgnet` object within a `mgnets`. The network data must be an `igraph` object
#' for `mgnet` objects. For `mgnets` objects, the network data should be a named list of `igraph` objects
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param value The new network data to be set.
#' @return The `mgnet` or `mgnets` object with updated network data.
#' @export
#' @importFrom methods validObject
#' @name netw<-
#' @aliases netw<-,mgnet-method netw<-,mgnets-method
setGeneric("netw<-", function(object, value) standardGeneric("netw<-"))

setMethod("netw<-", "mgnet", function(object, value) {
  
  # --- Ensure link_id exists if edges are present ----------------------------#
  value <- .ensure_link_id(value)
  
  object@netw <- value
  validObject(object)
  object
})

setMethod("netw<-", "mgnets", function(object, value) {
  is_list_mgnets_assign(object, value)
  for(i in names(object)) { netw(object@mgnets[[i]]) <- value[[i]] }
  validObject(object)
  object
})

# LINK (SET)
#------------------------------------------------------------------------------#
#' Update Network Links (Edges) for mgnet and mgnets Objects
#'
#' This function sets the network links (edges) for `mgnet` objects and also
#' updates each `mgnet` inside an `mgnets`. For `mgnet`, the value must be a
#' data frame / tibble with at least columns `from` and `to`. For `mgnets`,
#' the value can be either:
#' - a **named list** of data frames (one per contained `mgnet`, names matching `names(object)`), or
#' - a **single** data frame with a routing column `mgnet` (plus `from` and `to`)
#'   which is split and dispatched to each element.
#'
#' The existing vertex set and graph directedness are preserved from the
#' current network. Edge endpoints must be a subset of the existing vertex
#' names; otherwise an error is raised.
#'
#' @param object An `mgnet` or `mgnets` object to be updated.
#' @param value  The new edge table(s): a data frame/tibble for `mgnet`, or
#'   a named list of data frames / a single data frame with `mgnet` for `mgnets`.
#'
#' @return The updated `mgnet` or `mgnets` object with new network links.
#' @export
#' @name link<-
#' @aliases link<-,mgnet-method link<-,mgnets-method
setGeneric("link<-", function(object, value) standardGeneric("link<-"))

# -- helpers ------------------------------------------------------------------#
.is_edge_df <- function(x) {
  (is.data.frame(x) || inherits(x, "tbl_df")) &&
    all(c("from", "to") %in% colnames(x))
}

# Rebuild graph preserving vertices + directedness; enforce endpoints subset
.rebuild_graph_with_edges <- function(object, edge_df) {
  if (miss_netw(object)) {
    cli::cli_abort("No existing network available to derive vertices/directedness from.")
  }
  old_g  <- netw(object, selected = FALSE)
  vtab   <- igraph::as_data_frame(old_g, what = "vertices")
  is_dir <- igraph::is_directed(old_g)
  
  if (!("name" %in% colnames(vtab))) {
    cli::cli_abort("Vertex table must contain a {.val name} column.")
  }
  endpoints <- unique(c(edge_df$from, edge_df$to))
  unknown   <- setdiff(endpoints, vtab$name)
  if (length(unknown)) {
    cli::cli_abort(c(
      "x" = "Some edge endpoints are not present among existing vertex names.",
      "i" = "Unknown ids: {toString(utils::head(unknown, 10))}{if (length(unknown) > 10) ' ...' else ''}"
    ))
  }
  
  igraph::graph_from_data_frame(d = edge_df, directed = is_dir, vertices = vtab)
}


setMethod("link<-", c("mgnet", "ANY"), function(object, value) {
  
  objectName <- deparse(substitute(object))
  valueName  <- deparse(substitute(value))
  
  if (!.is_edge_df(value)) {
    cli::cli_abort("{.arg {valueName}} must be a tibble or data frame with columns {.val from} and {.val to}.")
  }
  if ("mgnet" %in% colnames(value)) {
    cli::cli_abort("{.arg {valueName}} cannot include a routing column {.val mgnet} when assigning to a single {.cls mgnet}.")
  }
  
  new_g <- .rebuild_graph_with_edges(object, value)
  object@netw <- new_g
  
  validObject(object)
  object
})


setMethod("link<-", c("mgnets","ANY"), function(object, value){
  
  if (is.data.frame(value) || inherits(value, "tbl_df")) {
    # single tibble/data.frame with `mgnet` routing column
    is_assign_mgnets_link_tbl(object, value)
    splitted_value <- split(value, value$mgnet)
    splitted_value <- lapply(splitted_value, \(x) { x$mgnet <- NULL; x })
    
    for (i in names(object)) link(object[[i]]) <- splitted_value[[i]]
    validObject(object)
    return(object)
    
  } else if (is.list(value)) {
    # named list path (like taxa<-)
    is_list_mgnets_assign(object, value)
    for (i in names(object)) link(object[[i]]) <- value[[i]]
    validObject(object)
    return(object)
    
  } else {
    cli::cli_abort("{.arg value} must be either a named list or a data.frame/tibble.")
  }
})



# COMMUNITY
#------------------------------------------------------------------------------#
#' Set Community Detection Results
#'
#' This setter function updates the community detection results for `mgnet` objects and
#' each `mgnet` object within a `mgnets`. The community data must be an object of class `communities`
#' for `mgnet` objects. For `mgnets` objects, the community data should be a named list of `communities` objects
#' corresponding to each `mgnet` object within the list.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param value The new community detection results to be set, which should be compatible
#' with the structure expected by the `mgnet` object's community slot. For individual `mgnet`
#' objects, this is typically an object of class `communities` as returned by community detection
#' functions in the `igraph` package. For `mgnets` objects, provide a named list where each element
#' is a `communities` object corresponding to the respective `mgnet` object within the list.
#'
#' @return The `mgnet` or `mgnets` object with updated community detection results.
#'
#' @export
#' @importFrom methods validObject
#' @name comm<-
#' @aliases comm<-,mgnet-method comm<-,mgnets-method
setGeneric("comm<-", function(object, value) standardGeneric("comm<-"))

setMethod("comm<-", "mgnet", function(object, value) {
  object@comm <- value
  validObject(object)
  object
})

setMethod("comm<-", "mgnets", function(object, value) {
  is_list_mgnets_assign(object, value)
  for(i in names(object)) { object@mgnets[[i]]@comm <- value[[i]] }
  validObject(object)
  object
})
