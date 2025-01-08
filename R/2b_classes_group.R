# TAXA
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Get Taxa Grouping Variables from an mgnet Object
#'
#' Retrieves the names of variables used for grouping taxa in an `mgnet` object.
#' These are set by [group_taxa()] and determine how subsequent functions
#' (e.g., `mutate_taxa`, `filter_taxa`, `mutate_netw`) operate on grouped taxa.
#'
#' @details
#' - If this function returns **`NULL`**, it indicates there is no explicit
#'   grouping information. In that scenario, expressions that do **not**
#'   reference abundance variables (`abun`, `rela`, `norm`) proceed **without** grouping;
#'   expressions that **do** reference those variables automatically group by
#'   `sample_id`.
#' - If this function returns an **empty character vector** (i.e., `character(0)`),
#'   it means grouping has been explicitly cleared with `group_taxa(object, NULL)`.
#'   No grouping is enforced, even for abundance-related expressions.
#' - If it returns a **character vector** of column names, those columns are
#'   used for grouping.
#' - For an `mgnetList`, each individual `mgnet` in the list may have its own
#'   grouping attributes. Typically, the `mgnet` identifier is always included
#'   in grouped operations across the list context.
#'
#' @param object An `mgnet` object or an `mgnetList`.
#'
#' @return 
#' - For a single `mgnet` object, a character vector of grouping variable names,
#'   or `NULL`/`character(0)` if no grouping is in effect.
#' - For an `mgnetList`, a list of such character vectors (one per `mgnet`).
#'
#' @export
setGeneric("get_group_taxa", function(object) standardGeneric("get_group_taxa"))

setMethod("get_group_taxa", "mgnet", function(object) {
  return(attr(object, "taxa_groups"))
})

setMethod("get_group_taxa", "mgnetList", function(object) {
  return(attr(object, "taxa_groups"))
})


#' Ungroup Taxa in an mgnet Object
#'
#' Clears the grouping attributes from the mgnet object, effectively removing any existing taxa grouping.
#' This action resets the context for operations in `mutate_taxa`, `filter_taxa`, and `mutate_netw`,
#' making them operate without any grouping constraints.
#'
#' @param object An `mgnet` object.
#'
#' @return The `mgnet` object without any grouping information.
#' @export
setGeneric("ungroup_taxa", function(object) standardGeneric("ungroup_taxa"))

setMethod("ungroup_taxa", "mgnet", function(object) {
  attr(object, "taxa_groups") <- NULL
  return(object)
})

setMethod("ungroup_taxa", "mgnetList", function(object) {
  attr(object, "taxa_groups") <- NULL
  return(object)
})

#' Group Taxa in an mgnet Object
#'
#' Groups taxa in an `mgnet` object by one or more specified variables from the taxa metadata.
#' Subsequent functions (e.g. `mutate_taxa`, `filter_taxa`, `mutate_netw`) will apply their operations
#' respecting these grouping variables.
#'
#' @details
#' - If you pass multiple variables (e.g., `Genus`, `Species`), taxa are grouped by the combination of those variables.
#' - If you pass exactly `NULL`, any existing taxa grouping is removed.
#' - **No Groups Specified**: If you do *not* specify any groups *and* your expressions do not reference abundance
#'   variables (`abun`, `rela`, `norm`), the transformations use **no grouping**. However, if your expressions **do**
#'   reference abundance variables, they automatically group by `sample_id` so that transformations apply per sample.
#' - **mgnetList Objects**: When grouping an `mgnetList`, the same logic applies to each `mgnet` inside the list.
#'
#' @param object An `mgnet` object (or `mgnetList`).
#' @param ... One or more unquoted expressions or column names from the taxa metadata.
#'   Each specified column must exist in the taxa data.  
#'   *To remove grouping, use `group_taxa(object, NULL)`.*
#'
#' @return The `mgnet` object (or `mgnetList`) with the new taxa grouping.
#' @importFrom rlang enquos as_name
#' @export
setGeneric("group_taxa", function(object, ...) standardGeneric("group_taxa"))

setMethod("group_taxa", "mgnet", function(object, ...) {
  group_vars <- rlang::enquos(...)
  
  # If exactly one argument is NULL, remove grouping
  if (length(group_vars) == 1 && rlang::quo_is_null(group_vars[[1]])) {
    attr(object, "taxa_groups") <- character(0)
    return(object)
  }
  
  # Otherwise, convert to character and validate
  var_names <- sapply(group_vars, rlang::as_name)
  missing_vars <- setdiff(var_names, taxa_vars(object))
  if (length(missing_vars)) {
    stop("Invalid meta variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Set grouping
  attr(object, "taxa_groups") <- var_names
  object
})

setMethod("group_taxa", "mgnetList", function(object, ...) {
  group_vars <- rlang::enquos(...)
  
  # If exactly one argument is NULL, remove grouping
  if (length(group_vars) == 1 && rlang::quo_is_null(group_vars[[1]])) {
    attr(object, "taxa_groups") <- character(0)
    return(object)
  }
  
  # Otherwise, convert to character and validate
  var_names <- sapply(group_vars, rlang::as_name)
  missing_vars <- setdiff(var_names, taxa_vars(object, "unique"))
  if (length(missing_vars)) {
    stop("Invalid meta variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Store the variable names as an attribute
  attr(object[[i]], "taxa_groups") <- var_names
  return(object)
})


# META
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Get Meta Grouping Variables from an mgnet Object
#'
#' Retrieves the names of variables used for grouping the meta data in an `mgnet` object.
#' These are set by [group_meta()] and determine how subsequent functions
#' (e.g., `mutate_meta`, `filter_meta`) operate on grouped meta.
#'
#' @details
#' - If this function returns **`NULL`**, there is no explicit grouping information.
#'   In that case, expressions that do **not** reference abundance variables
#'   (`abun`, `rela`, `norm`) have **no grouping**, while expressions referencing them
#'   default to grouping by `sample_id`.
#' - If it returns an **empty character vector** (i.e., `character(0)`),
#'   grouping was explicitly removed via `group_meta(object, NULL)`. No grouping is
#'   enforced, even for abundance-related expressions.
#' - If it returns a **character vector** of column names, those columns
#'   are used for grouping meta data.
#' - For an `mgnetList`, each `mgnet` may store its own grouping variables, and
#'   combined operations typically account for each `mgnet` in the list.
#'
#' @param object An `mgnet` object or an `mgnetList`.
#'
#' @return 
#' - For a single `mgnet`, a character vector of grouping variable names, or
#'   `NULL`/`character(0)` if no grouping is in effect.
#' - For an `mgnetList`, a list of such character vectors (one per `mgnet`).
#'
#' @export
setGeneric("get_group_meta", function(object) standardGeneric("get_group_meta"))

setMethod("get_group_meta", "mgnet", function(object) {
  return(attr(object, "meta_groups"))
})

setMethod("get_group_meta", "mgnetList", function(object) {
  return(attr(object, "meta_groups"))
})


#' Ungroup Meta in an mgnet Object
#'
#' Clears the grouping attributes from the mgnet object, effectively removing any existing meta grouping.
#' This action resets the context for operations in `mutate_meta`, `filter_meta`,
#' making them operate without any grouping constraints.
#'
#' @param object An `mgnet` object.
#'
#' @return The `mgnet` object without any grouping information.
#' @export
setGeneric("ungroup_meta", function(object) standardGeneric("ungroup_meta"))

setMethod("ungroup_meta", "mgnet", function(object) {
  attr(object, "meta_groups") <- NULL
  return(object)
})

setMethod("ungroup_meta", "mgnetList", function(object) {
  attr(object, "meta_groups") <- NULL
  return(object)
})


#' Group Meta in an mgnet Object
#'
#' Groups the meta data in an `mgnet` object by one or more specified variables from the meta data.
#' Functions like `mutate_meta` and `filter_meta` will apply their operations based on these grouping variables.
#'
#' @details
#' - If you pass multiple variables (e.g., `Site`, `Treatment`), the meta data are grouped by the combination of those variables.
#' - If you pass exactly `NULL`, any existing meta grouping is removed.
#' - **No Groups Specified**: If you do *not* specify any groups *and* your expressions do not reference abundance
#'   variables (`abun`, `rela`, `norm`), the transformations use **no grouping**. However, if your expressions **do**
#'   reference abundance variables, they automatically group by `sample_id` so that transformations apply per sample.
#' - **mgnetList Objects**: When grouping an `mgnetList`, the same logic applies to each `mgnet` inside the list.
#'
#' @param object An `mgnet` object (or `mgnetList`).
#' @param ... One or more unquoted expressions or column names from the meta data.
#'   Each specified column must exist in the meta data.  
#'   *To remove grouping, use `group_meta(object, NULL)`.*
#'
#' @return The `mgnet` object (or `mgnetList`) with the new meta grouping.
#' @importFrom rlang enquos as_name
#' @export
setGeneric("group_meta", function(object, ...) standardGeneric("group_meta"))

setMethod("group_meta", "mgnet", function(object, ...) {
  group_vars <- rlang::enquos(...)
  
  # If exactly one argument is NULL, remove grouping
  if (length(group_vars) == 1 && rlang::quo_is_null(group_vars[[1]])) {
    attr(object, "meta_groups") <- character(0)
    return(object)
  }
  
  # Otherwise, convert to character and validate
  var_names <- sapply(group_vars, rlang::as_name)
  missing_vars <- setdiff(var_names, meta_vars(object))
  if (length(missing_vars)) {
    stop("Invalid meta variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Set grouping
  attr(object, "meta_groups") <- var_names
  object
})

setMethod("group_meta", "mgnetList", function(object, ...) {
  group_vars <- rlang::enquos(...)
  
  # If exactly one argument is NULL, remove grouping
  if (length(group_vars) == 1 && rlang::quo_is_null(group_vars[[1]])) {
    attr(object, "meta_groups") <- character(0)
    return(object)
  }
  
  # Otherwise, convert to character and validate
  var_names <- sapply(group_vars, rlang::as_name)
  missing_vars <- setdiff(var_names, meta_vars(object, "unique"))
  if (length(missing_vars)) {
    stop("Invalid meta variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Store the variable names as an attribute
  attr(object, "meta_groups") <- var_names
  return(object)
})
