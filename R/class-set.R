#' Set slots of an `mgnet` object using tidy-style expressions
#'
#' @description
#' `set_mgnet()` lets you compute and assign one or more slots of an `mgnet`
#' object using rlang-style expressions that can reference existing slots by
#' their bare names (e.g., `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`).
#' Expressions are evaluated **sequentially** in the order provided, so later
#' expressions can depend on values assigned by earlier ones.
#'
#' @param object An object of class `mgnet`.
#' @param ... Named expressions, where each name must be one of the allowed
#'   slot names: `abun`, `rela`, `norm`, `meta`, `taxa`, `netw`, `comm`.
#' @param .cores Integer number of CPU cores to use when `object` is an
#'   `mgnets` object. If not supplied, the method for `mgnets` will use
#'   `min(length(object), parallel::detectCores() - 1)` (with a safe
#'   fallback to 1). For `mgnet` objects, `.cores` is ignored and computation
#'   is always sequential.
#' @param .packages Optional character vector of package names to be loaded
#'   on worker nodes when `.cores > 1`. This is useful when expressions
#'   in `set_mgnet()` use functions not qualified with `pkg::fun`. Packages
#'   must already be installed. Defaults include `"mgnet"` and `"magrittr"`.
#'
#' @return The modified `mgnet` object with updated slots.
#'
#' @section How name resolution works:
#' Within each expression, the symbols `abun`, `rela`, `norm`, `meta`, `taxa`,
#' `netw`, and `comm` resolve to the **current** working values (initially taken
#' from the object, then updated as each expression is assigned).
#'
#' @export
setGeneric("set_mgnet", function(object, ..., .cores, .packages = NULL) standardGeneric("set_mgnet"))

#' @rdname set_mgnet
#' @export
setMethod("set_mgnet", signature(object = "mgnet"), function(object, ..., .cores, .packages = NULL) {
  # Capture user expressions as quosures, preserving their environment
  dots <- rlang::enquos(..., .named = TRUE)
  
  # Allowed target slots
  allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
  
  # 1) Basic validations on `...`
  if (length(dots) == 0L) {
    cli::cli_abort("No expressions supplied to {.fn set_mgnet}.")
  }
  
  nm <- names(dots)
  if (is.null(nm) || any(nm == "")) {
    cli::cli_abort("All arguments to {.fn set_mgnet} must be named.")
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown)) {
    cli::cli_abort(
      c(
        "Some targets are not valid {'.mgnet'} slots:",
        "x" = paste0(unknown, collapse = ", "),
        "i" = "Valid targets are: {paste(allowed, collapse = ', ')}"
      )
    )
  }
  
  # 2) Build the initial data mask from current object values
  mask <- list(
    abun = abun(object),
    rela = rela(object),
    norm = norm(object),
    meta = meta(object),
    taxa = taxa(object),
    netw = netw(object),
    comm = comm(object)
  )
  
  # 3) Evaluate each expression with rlang::eval_tidy in the mask
  for (target in nm) {
    quo <- dots[[target]]
    
    value <- tryCatch(
      rlang::eval_tidy(quo, data = mask, env = rlang::get_env(quo)),
      error = function(e) {
        cli::cli_abort(c(
          "Failed while computing {.field {target}} in {.fn set_mgnet}.",
          "x" = conditionMessage(e),
          "i" = "Within expressions you can reference: {paste(allowed, collapse = ', ')}"
        ))
      }
    )
    
    # Assign via your replacement accessors and refresh the mask
    object <- do.call(paste0(target, "<-"), list(object, value))
    mask[[target]] <- value
  }
  
  # 4) Return the updated object (validity is handled inside your setters)
  object
})


#' @rdname set_mgnet
#' @export
setMethod(
  "set_mgnet",
  signature(object = "mgnets"),
  function(object, ..., .cores, .packages = NULL) {
    
    dots <- rlang::enquos(..., .named = TRUE)
    allowed <- c("abun", "rela", "norm", "meta", "taxa", "netw", "comm")
    
    if (length(dots) == 0L) {
      cli::cli_abort("No expressions supplied to {.fn set_mgnet}.")
    }
    
    nm <- names(dots)
    if (is.null(nm) || any(nm == "")) {
      cli::cli_abort("All arguments to {.fn set_mgnet} must be named.")
    }
    
    unknown <- setdiff(nm, allowed)
    if (length(unknown)) {
      cli::cli_abort(
        c(
          "Some targets are not valid {'.mgnet'} slots:",
          "x" = paste0(unknown, collapse = ", "),
          "i" = "Valid targets are: {paste(allowed, collapse = ', ')}"
        )
      )
    }
    
    n <- length(object)
    if (n == 0L) return(object)
    
    ## ---- .cores validation ----
    if (missing(.cores) || is.null(.cores)) {
      # not specified: use min(n, detectCores() - 1)
      cores_sys <- parallel::detectCores()
      if (is.na(cores_sys) || cores_sys <= 1L) {
        .cores <- 1L
      } else {
        .cores <- min(n, cores_sys - 1L)
      }
    } else {
      # user-specified: must be a single integer >= 1
      if (!is.numeric(.cores) || length(.cores) != 1L || is.na(.cores)) {
        cli::cli_abort(
          "{.arg .cores} must be a single, non-missing numeric value."
        )
      }
      if (.cores < 1) {
        cli::cli_abort(
          "{.arg .cores} must be an integer greater or equal to 1."
        )
      }
      if (.cores %% 1 != 0) {
        cli::cli_abort(
          "{.arg .cores} must be an integer (no decimal part)."
        )
      }
      # do not use more workers than elements
      if (.cores > n) {
        .cores <- n
      }
    }
    
    .cores <- as.integer(.cores)
    
    ## Helper that runs set_mgnet() on a single element
    worker <- function(i) {
      el_name <- names(object)[i]
      label <- if (!is.null(el_name) && nzchar(el_name)) {
        paste0("element '", el_name, "' (index ", i, ")")
      } else {
        paste0("element at index ", i)
      }
      
      om <- object[[i]]
      if (!methods::is(om, "mgnet")) {
        cli::cli_abort(c(
          "All elements of an {.cls mgnets} object must be {.cls mgnet}.",
          "x" = paste0(
            "Found element of class {.cls ", class(om)[1],
            "} at index ", i, "."
          )
        ))
      }
      
      # run set_mgnet() on the single element
      tryCatch(
        rlang::inject(set_mgnet(om, !!!dots)),
        error = function(e) {
          cli::cli_abort(c(
            "Failed while running {.fn set_mgnet} on {.val {label}}.",
            "x" = conditionMessage(e)
          ))
        }
      )
    }
    
    idx <- seq_len(n)
    
    if (.cores == 1L) {
      ## --- sequential version (original behavior) ---
      for (i in idx) {
        object[[i]] <- worker(i)
      }
      return(object)
    }
    
    ## --- parallel version with PSOCK cluster (Windows-friendly) ---
    cl <- parallel::makeCluster(.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Build list of packages to load on workers
    pkgs_to_load <- unique(c("mgnet", "magrittr", .packages))
    
    # Validate packages before sending to workers
    for (p in pkgs_to_load) {
      if (!requireNamespace(p, quietly = TRUE)) {
        cli::cli_abort("Package {.val {p}} in {.arg .packages} is not installed.")
      }
    }
    
    # Export pkgs_to_load to workers so it is visible inside clusterEvalQ()
    parallel::clusterExport(
      cl,
      varlist = "pkgs_to_load",
      envir = environment()
    )
    
    # Load packages in each worker
    parallel::clusterEvalQ(cl, {
      for (p in pkgs_to_load) {
        suppressPackageStartupMessages(
          require(p, character.only = TRUE)
        )
      }
      NULL
    })
    
    # parLapply serializes 'worker', 'object' and 'dots' and sends them to workers
    res <- parallel::parLapply(cl, idx, worker)
    
    for (i in idx) {
      object[[i]] <- res[[i]]
    }
    
    object
  }
)
