#' @include class-mgnet.R class-mgnets.R
NULL

#==============================================================================#
# Merge mgnet / mgnets
#==============================================================================#

#' Merge `mgnet` or `mgnets` objects
#'
#' @description
#' Merge multiple `mgnet` objects into a single `mgnet`, or merge multiple
#' `mgnets` objects name-by-name into a single `mgnets`.
#'
#' Input objects must be all of class `mgnet` or all of class `mgnets`.
#' Mixed input is not allowed.
#'
#' For each slot, the user may supply a merge function. If a merge function is
#' not provided for a given slot, that slot is left empty in the merged result.
#'
#' @param ... `mgnet` objects, `mgnets` objects, or lists of such objects.
#' @param abun_fun,rela_fun,norm_fun Functions used to merge abundance-like slots.
#' @param meta_fun,taxa_fun Functions used to merge metadata tables.
#' @param netw_fun,comm_fun Functions used to merge network and community slots.
#'
#' @return
#' An object of class `mgnet` if the inputs are `mgnet` objects, or an object
#' of class `mgnets` if the inputs are `mgnets` objects.
#'
#' @name merge_mgnet
#' @export
merge_mgnet <- function(...,
                        abun_fun = NULL,
                        rela_fun = NULL,
                        norm_fun = NULL,
                        meta_fun = NULL,
                        taxa_fun = NULL,
                        netw_fun = NULL,
                        comm_fun = NULL) {
  
  check_fun <- function(x, name) {
    if (!is.null(x) && !is.function(x)) {
      cli::cli_abort("{.arg {name}} must be a function.")
    }
  }
  
  check_fun(abun_fun, "abun_fun")
  check_fun(rela_fun, "rela_fun")
  check_fun(norm_fun, "norm_fun")
  check_fun(meta_fun, "meta_fun")
  check_fun(taxa_fun, "taxa_fun")
  check_fun(netw_fun, "netw_fun")
  check_fun(comm_fun, "comm_fun")
  
  objs <- list(...)
  
  bad_named <- names(objs)
  bad_named <- bad_named[!is.null(bad_named) & nzchar(bad_named)]
  
  if (length(bad_named) > 0L) {
    cli::cli_abort(c(
      "x" = "Unexpected named arguments were supplied in `{...}`.",
      "i" = "Merge functions must be passed using the `*_fun` arguments.",
      "v" = "Unexpected names: {.val {unique(bad_named)}}."
    ))
  }
  
  objs <- unlist(objs, recursive = FALSE)
  
  if (length(objs) == 0L) {
    cli::cli_abort("No objects supplied.")
  }
  
  are_mgnet  <- all(vapply(objs, methods::is, logical(1), class2 = "mgnet"))
  are_mgnets <- all(vapply(objs, methods::is, logical(1), class2 = "mgnets"))
  
  if (!(are_mgnet || are_mgnets)) {
    cli::cli_abort(
      c(
        "x" = "All inputs must be either {.cls mgnet} objects or {.cls mgnets} objects.",
        "i" = "Mixed input is not supported."
      )
    )
  }
  
  if (length(objs) == 1L) {
    return(objs[[1L]])
  }
  
  #--------------------------------------------------------------------------#
  # Merge mgnet objects
  #--------------------------------------------------------------------------#
  if (are_mgnet) {
    
    abun_merged <- matrix(numeric(0), nrow = 0, ncol = 0)
    rela_merged <- matrix(numeric(0), nrow = 0, ncol = 0)
    norm_merged <- matrix(numeric(0), nrow = 0, ncol = 0)
    meta_merged <- data.frame()
    taxa_merged <- data.frame()
    netw_merged <- igraph::make_empty_graph(0, directed = FALSE)
    comm_merged <- NULL
    
    if (!is.null(abun_fun)) {
      abun_merged <- purrr::map(objs, \(x) abun(x, .fmt = "df")) %>%
        purrr::reduce(abun_fun) %>%
        as.matrix()
    }
    
    if (!is.null(rela_fun)) {
      rela_merged <- purrr::map(objs, \(x) rela(x, .fmt = "df")) %>%
        purrr::reduce(rela_fun) %>%
        as.matrix()
    }
    
    if (!is.null(norm_fun)) {
      norm_merged <- purrr::map(objs, \(x) norm(x, .fmt = "df")) %>%
        purrr::reduce(norm_fun) %>%
        as.matrix()
    }
    
    if (!is.null(meta_fun)) {
      meta_merged <- purrr::map(objs, \(x) meta(x, .fmt = "df")) %>%
        purrr::reduce(meta_fun)
    }
    
    if (!is.null(taxa_fun)) {
      taxa_merged <- purrr::map(objs, \(x) taxa(x, .fmt = "df")) %>%
        purrr::reduce(taxa_fun)
    }
    
    if (!is.null(netw_fun)) {
      netw_merged <- purrr::map(objs, \(x) netw(x, selected = FALSE)) %>%
        purrr::reduce(netw_fun)
    }
    
    if (!is.null(comm_fun)) {
      comm_merged <- purrr::map(objs, \(x) comm(x)) %>%
        purrr::reduce(comm_fun)
    }
    
    return(
      mgnet(
        abun = abun_merged,
        rela = rela_merged,
        norm = norm_merged,
        meta = meta_merged,
        taxa = taxa_merged,
        netw = netw_merged,
        comm = comm_merged
      )
    )
  }
  
  #--------------------------------------------------------------------------#
  # Merge mgnets objects
  #--------------------------------------------------------------------------#
  names_list <- purrr::map(objs, names)
  are_names_equal <- purrr::reduce(names_list, \(a, b) setequal(a, b))
  
  if (!are_names_equal) {
    cli::cli_abort(
      "All {.cls mgnets} objects must contain the same set of names."
    )
  }
  
  merged_names <- names(objs[[1]])
  results <- stats::setNames(vector("list", length(merged_names)), merged_names)
  
  for (nm in merged_names) {
    results[[nm]] <- do.call(
      merge_mgnet,
      c(
        lapply(objs, `[[`, nm),
        list(
          abun_fun = abun_fun,
          rela_fun = rela_fun,
          norm_fun = norm_fun,
          meta_fun = meta_fun,
          taxa_fun = taxa_fun,
          netw_fun = netw_fun,
          comm_fun = comm_fun
        )
      )
    )
  }
  
  mgnets(results)
}


#==============================================================================#
# Collapse groups of mgnet objects inside an mgnets
#==============================================================================#

#' Collapse groups of `mgnet` objects within an `mgnets`
#'
#' @description
#' Collapse subsets of an `mgnets` object into merged `mgnet` objects according
#' to a user-specified grouping.
#'
#' @param object An object of class `mgnets`.
#' @param by A named or unnamed list defining groups of `mgnet` names to merge.
#' Each element must be a character vector of names found in `object`.
#' @param abun,rela,norm,meta,taxa,netw,comm Functions used to merge the
#' corresponding slots.
#'
#' @return
#' An object of class `mgnets`, where each output element is the merge of one
#' group specified in `by`.
#'
#' @details
#' If `by` is unnamed, names are generated by concatenating the grouped element
#' names with `"-"`.
#'
#' Output group names must be unique.
#'
#' @name mgnet_collapse
#' @export
mgnet_collapse <- function(object,
                           by,
                           abun = NULL,
                           rela = NULL,
                           norm = NULL,
                           meta = NULL,
                           taxa = NULL,
                           netw = NULL,
                           comm = NULL) {
  
  if (!methods::is(object, "mgnets")) {
    cli::cli_abort("{.arg object} must be an object of class {.cls mgnets}.")
  }
  
  if (!is.list(by)) {
    cli::cli_abort("{.arg by} must be a list of character vectors.")
  }
  
  if (length(by) == 0L) {
    cli::cli_abort("{.arg by} must contain at least one group.")
  }
  
  if (!all(vapply(by, is.character, logical(1)))) {
    cli::cli_abort("Each element of {.arg by} must be a character vector.")
  }
  
  by_names <- names(by)
  if (is.null(by_names)) {
    by_names <- rep("", length(by))
  }
  
  names(by) <- ifelse(
    nzchar(by_names),
    by_names,
    vapply(by, function(x) paste(x, collapse = "-"), character(1))
  )
  
  if (anyDuplicated(names(by)) > 0L) {
    cli::cli_abort("Output group names in {.arg by} must be unique.")
  }
  
  object_names <- names(object)
  
  unknown <- unique(setdiff(unlist(by, use.names = FALSE), object_names))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Some names in {.arg by} are not present in {.arg object}.",
        "v" = "Unknown names: {.val {unknown}}."
      )
    )
  }
  
  results <- stats::setNames(vector("list", length(by)), names(by))
  object_list <- as.list(object)
  
  for (i in seq_along(by)) {
    sub_object <- object_list[by[[i]]]
    
    results[[i]] <- do.call(
      merge_mgnet,
      c(
        sub_object,
        list(
          abun = abun,
          rela = rela,
          norm = norm,
          meta = meta,
          taxa = taxa,
          netw = netw,
          comm = comm
        )
      )
    )
  }
  
  mgnets(results)
}