# SPLIT MGNET
#------------------------------------------------------------------------------#
#' Split an mgnet Object into an mgnetList Based on Specified Columns
#'
#' @description
#' This method splits an `mgnet` object into subsets based on unique combinations
#' of values in specified columns from the `info_sample` data frame. Each subset
#' corresponds to an `mgnet` object in the resulting `mgnetList`, facilitating
#' analyses specific to each unique combination of metadata values.
#'
#' @param object An `mgnet` object.
#' @param ... Columns in the `info_sample` data frame used to define subsets.
#'        You can specify columns by name, using any of the selection methods supported by `dplyr::select()`.
#' @return An `mgnetList` object containing `mgnet` objects for each unique combination
#'         of specified column values.
#'
#' @details
#' The `split_mgnet` function leverages `dplyr` functionality to identify unique combinations
#' of specified column values in the `info_sample` data frame. It then creates a new `mgnet`
#' object for each combination, allowing for tailored analysis per group.
#'
#' This function is particularly useful in exploratory data analysis and preprocessing stages
#' where data need to be examined or analyzed based on specific grouping variables.
#'
#' @note
#' This function depends on `dplyr` for data manipulation. Ensure that `dplyr` is installed
#' and loaded in your R session.
#'
#' @export
#' @aliases split_mgnet,mgnet-method
#' @importFrom dplyr select distinct group_by group_split semi_join
#' @importFrom rlang enquos
setGeneric("split_mgnet", function(object, ...) {
  standardGeneric("split_mgnet")
})

setMethod("split_mgnet", "mgnet", function(object, ...) {
  # Capture the column names as quosures
  split_cols <- enquos(...)
  
  # Check if split_cols is empty
  if (length(split_cols) == 0) {
    stop("No columns specified for splitting. Please specify one or more columns.")
  }

  # Get the unique combinations from info_sample for the specified columns
  info_sample <- info_sample(object, "tbl")
  groups <- info_sample %>%
    dplyr::select(!!!split_cols) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!!split_cols) %>%
    dplyr::group_split()

  # Initialize an empty list to store the subsetted mgnet objects
  mgnet_list <- mgnetList()

  # For each group, filter the info_sample based on the unique combination and create a new mgnet object
  for (i in seq_along(groups)) {
    # Filter info_sample based on the group
    group_df <- groups[[i]]
    filtered_info_sample <- dplyr::semi_join(info_sample, group_df, by = names(group_df))

    # Filter the original mgnet object to create a subset based on the filtered info_sample
    subsetted_mgnet <- object[filtered_info_sample$sample_id, , drop = FALSE]

    # Add the subsetted mgnet object to the list with a meaningful name
    group_name <- apply(group_df, 1, function(row) {paste(row, collapse = "-")})
    mgnet_list[[group_name]] <- subsetted_mgnet
  }

  return(mgnet_list)
})