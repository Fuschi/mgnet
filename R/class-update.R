#' @include class-mgnet.R class-mgnets.R
NULL

#==============================================================================#
# Update slots of `mgnet` / `mgnets`
#==============================================================================#

#' Update slots of an `mgnet` or `mgnets` object
#'
#' @description
#' Evaluate one or more named expressions sequentially and use their results to
#' update slots of an `mgnet` object. Supported targets are:
#' `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, and `comm`.
#'
#' Each expression is evaluated in a mask containing the current values of these
#' same targets. Updates are applied in order, so later expressions can depend
#' on earlier ones.
#'
#' For `mgnets`, the same update is applied to each contained `mgnet` object,
#' optionally in parallel.
#'
#' @param object An `mgnet` or `mgnets` object.
#' @param ... Named expressions. The names indicate which target slot to update.
#' @param .cores Number of worker processes used for `mgnets`. Defaults to `1`.
#' @param .packages Optional character vector of package names to load on worker
#'   processes when `.cores > 1`.
#'
#' @details
#' Expressions are evaluated sequentially. After each successful assignment, the
#' corresponding value is re-read from the updated object and inserted back into
#' the evaluation mask. This ensures that subsequent expressions see the actual
#' stored value after any setter-side reordering, coercion, or normalization.
#'
#' Supported targets are:
#' \itemize{
#'   \item `abun`
#'   \item `rela`
#'   \item `norm`
#'   \item `meta`
#'   \item `taxa`
#'   \item `netw`
#'   \item `comm`
#' }
#'
#' @return An updated `mgnet` or `mgnets` object.
#'
#' @name update_mgnet
NULL

#' @rdname update_mgnet
#' @export
setGeneric("update_mgnet", function(object, ..., .cores = 1L, .packages = NULL) {
  standardGeneric("update_mgnet")
})

#------------------------------------------------------------------------------#
# mgnet method
#------------------------------------------------------------------------------#

#' @rdname update_mgnet
#' @export
setMethod("update_mgnet", "mgnet", function(object, ..., .cores = 1L, .packages = NULL) {
  
  dots <- rlang::enquos(..., .named = TRUE)
  
  if (length(dots) == 0L) {
    return(object)
  }
  
  nm <- names(dots)
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  unknown <- setdiff(nm, allowed)
  
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Some targets are not valid slots of {.cls mgnet}: {.val {unknown}}.",
        "i" = "Valid targets are: {.val {allowed}}."
      )
    )
  }
  
  # Initial evaluation mask
  mask <- list(
    abun = abun(object),
    rela = rela(object),
    norm = norm(object),
    meta = meta(object),
    taxa = taxa(object),
    netw = netw(object),
    comm = comm(object)
  )
  
  for (target in nm) {
    expr <- dots[[target]]
    
    value <- rlang::eval_tidy(expr, data = mask, env = rlang::caller_env())
    
    object <- do.call(paste0(target, "<-"), list(object, value))
    
    # Refresh mask with the actual stored value after setter-side processing
    mask[[target]] <- switch(
      target,
      abun = abun(object),
      rela = rela(object),
      norm = norm(object),
      meta = meta(object),
      taxa = taxa(object),
      netw = netw(object),
      comm = comm(object)
    )
  }
  
  methods::validObject(object)
  object
})

#------------------------------------------------------------------------------#
# mgnets method
#------------------------------------------------------------------------------#

#' @rdname update_mgnet
#' @export
setMethod("update_mgnet", "mgnets", function(object, ..., .cores = 1L, .packages = NULL) {
  
  dots <- rlang::enquos(..., .named = TRUE)
  
  if (length(dots) == 0L) {
    return(object)
  }
  
  n <- length(object)
  if (n == 0L) {
    return(object)
  }
  
  if (!is.numeric(.cores) || length(.cores) != 1L || is.na(.cores) || .cores < 1) {
    cli::cli_abort("{.arg .cores} must be a single positive integer.")
  }
  
  .cores <- as.integer(.cores)
  .cores <- min(.cores, n)
  
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  nm <- names(dots)
  unknown <- setdiff(nm, allowed)
  
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Some targets are not valid slots of {.cls mgnet}: {.val {unknown}}.",
        "i" = "Valid targets are: {.val {allowed}}."
      )
    )
  }
  
  if (.cores == 1L) {
    for (i in names(object)) {
      om <- object[[i]]
      
      if (!methods::is(om, "mgnet")) {
        cli::cli_abort("Element {.val {i}} is not an {.cls mgnet} object.")
      }
      
      object[[i]] <- rlang::inject(update_mgnet(om, !!!dots))
    }
    
    methods::validObject(object)
    return(object)
  }
  
  cl <- parallel::makePSOCKcluster(.cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (!is.null(.packages)) {
    parallel::clusterEvalQ(cl, {
      invisible(lapply(.packages, require, character.only = TRUE))
    })
  }
  
  elements <- as.list(object)
  
  worker <- function(om, dots) {
    if (!methods::is(om, "mgnet")) {
      stop("Worker received a non-mgnet object.")
    }
    rlang::inject(update_mgnet(om, !!!dots))
  }
  
  res <- parallel::parLapply(
    cl,
    X = elements,
    fun = worker,
    dots = dots
  )
  
  names(res) <- names(object)
  object@mgnets <- res
  
  methods::validObject(object)
  object
})