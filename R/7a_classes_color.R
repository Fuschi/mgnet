# SET CATEGORICAL COLOR
#------------------------------------------------------------------------------#
#' Set Categorical Colors in `mgnet` or `mgnetList` Objects
#'
#' @description
#' This function assigns colors to taxa based on the results of specified expressions evaluated within the context
#' of defined grouping variables. It allows for the dynamic creation of color-coded categories for visual representation
#' in subsequent visualizations.
#'
#' @param object An `mgnet` or `mgnetList` object containing the taxa data.
#' @param ... The expressions to evaluate, which determine how colors are assigned within each grouping variable's context.
#'        Each expression must be named as it dictates the resulting column name in the taxa metadata.
#' @param vars Character vector specifying one or more grouping variables by which taxa are categorized before coloring.
#'        Each variable specified must be present in the taxa data.
#' @param decreasing Logical; whether to sort the categories within each group in decreasing order based on the expression evaluation.
#' @param distinct_colors Integer; the maximum number of distinct colors to generate for the categories. Defaults to 20.
#' @param colorspace Character or list; the color palette to use, which can be a predefined palette name or a list specifying
#'        hue, saturation, and lightness ranges. Supported names are 'pretty', 'pretty_dark', 'rainbow', and 'pastels'.
#'        For a detailed explanation of these options, see the `qualpalr` package documentation or enter `??qualpalr` in the R console.
#' @param alpha Numeric; transparency level of the colors, where 1 is opaque and 0 is fully transparent. Defaults to 1.
#' @param extraColor Character; the color used for any additional categories beyond `distinct_colors`. Defaults to white ("#FFFFFF").
#'
#' @details
#' `set_categorical_color` dynamically generates colors based on the results of specified expressions, applying these colors
#' to the taxa metadata. This method is useful for preparing data for visualization, where colors can help in distinguishing
#' taxa based on categorical attributes derived from the data.
#'
#' The function modifies the original `mgnet` or `mgnetList` object by adding new columns that represent color codes,
#' corresponding to each grouping variable and expression combination.
#'
#' If any existing columns in the taxa data conflict with the newly generated color column names, they will be overwritten
#' with the new color data, and a warning will be issued.
#'
#' @return The modified `mgnet` or `mgnetList` object with added color columns in the taxa data.
#'
#' @seealso
#' \code{\link[qualpalr]{qualpalr}} for details on generating qualitative color palettes.
#'
#' @examples
#' # Assuming `mg` is an `mgnet` object with taxa metadata containing `genus` and abundance data
#' mg <- set_categorical_color(mg, sum_abun = sum(abun), vars = "genus")
#' @export
#' @importFrom dplyr mutate group_by arrange
#' @importFrom tidyr expand_grid
#' @importFrom purrr map_chr
#' @importFrom methods slot
#' @importFrom tibble as_tibble
#' @export
#' @importFrom dplyr mutate %>%
#' @aliases set_categorical_color,mgnet-method set_categorical_color,mgnetList-method
setGeneric("set_categorical_color", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                             distinct_colors = 20, colorspace = "pretty",
                                             alpha = 1, extraColor = "#FFFFFF"
                                             ) standardGeneric("set_categorical_color"))

setMethod("set_categorical_color", signature = "mgnet", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                                                 distinct_colors = 20, colorspace = "pretty",
                                                                 alpha = 1, extraColor = "#FFFFFF"){

    # CHECKS
    #----------------------------------------------------------------------------#
  
    # Capture all the expressions provided
    expressions <- rlang::enquos(...)
    expr_names <- names(expressions)
    
    # Check if any of the expressions are unnamed
    if (any(nzchar(expr_names) == 0)) {
      stop("All the expressions must be named.")
    }

    # Check the reserved keywords
    check_reserved_keywords(expressions)

    # Check all expressions are named
    if(any(nzchar(expr_names)==0) | any(duplicated(expr_names))){
      stop("All the expressions must have an unique name")
    }

    # Capture required keys from expressions
    keys_required <- expressions %>%
      purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
      unlist() %>%
      unique()

    # Store needed abundances keys
    needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))

    # Check the variables needed
    validate_required_variables(object, keys_required, "taxa")

    if (!is.character(vars) || "sample_id" %in% vars) {
      stop("Error: 'vars' must be a character vector and cannot include 'sample_id'.")
    }

    # Forbidden functions and disallowed variables
    check_forbidden_expressions(expressions)
    
    # Check the variables needed
    validate_required_variables(object, keys_required, "taxa")
    
    # Generate new column names by combining variable elements with expression names
    new_column_names <- purrr::map2_chr(vars, expr_names, ~ paste("color", .x, .y, sep = "_"))
    
    # Check if these new column names already exist in the taxa
    existing_columns <- colnames(taxa(object, .fmt = "tbl"))
    if (any(new_column_names %in% existing_columns)) {
      conflicting_names <- new_column_names[new_column_names %in% existing_columns]
      warning("The following column names in taxa data will be overwritten: ", paste(conflicting_names, collapse = ", "))
    }
    
    if(distinct_colors > 99 || distinct_colors < 0){
      stop("Error: Ddstinct_colors must be in range [0,99].")
    }
    
    if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
      stop("Error: alpha must be in the range [0,1].")
    }
    
    if(!is.character(colorspace) && !is.list(colorspace)) {
      stop("Error: colorspace must be a character or a list.")
    }
    
    if(is.character(colorspace)){
      if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
        stop("Error: when colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
      }
    }
    
    if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
      stop("Error: when colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpalr.")
    }
    
    if(!is.logical(decreasing)) stop("Error: decreasing must be logical.")
    
    # Compile list of valid variables including "taxa_id" and additional taxa variables
    valid_vars <- c("taxa_id", taxa_vars(object))
    
    # Check if all elements of 'var' are in 'valid_vars'
    if (!all(vars %in% valid_vars)) {
      missing_vars <- vars[!vars %in% valid_vars]
      stop("Error: The following grouping 'vars' elements are not valid or missing in the taxa data: ", 
           paste(missing_vars, collapse = ", "))
    }

    # END CHECKS
    #----------------------------------------------------------------------------#

    # CREATE THE BASE FOR THE SOLUTION
    #----------------------------------------------------------------------------#
    # Initialize sample_info_mutated
    if (length(taxa(object)) != 0) {
      taxa_info_mutated <- taxa(object, .fmt = "tbl")
    } else {
      taxa_info_mutated <- tibble::tibble(taxa_id = taxa_id(object))
    }

    # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
    #----------------------------------------------------------------------------#
    if (length(needed_abundance_keys) > 0) {
      # Create the grid of sample_id and taxa_id combinations
      long_abun <- tidyr::expand_grid(
        sample_id = sample_id(object),
        taxa_id = taxa_id(object)
      )

      # Left join all abundance variables needed
      for (abundance_key in needed_abundance_keys) {
        long_abun <- long_abun %>%
          dplyr::left_join(
            methods::slot(object, abundance_key) %>%
              tibble::as_tibble(rownames = "sample_id") %>%
              tidyr::pivot_longer(-sample_id,
                                  names_to = "taxa_id",
                                  values_to = abundance_key),
            by = c("sample_id", "taxa_id")
          )
      }
      
      if(length(taxa(object))!=0) long_abun <- dplyr::left_join(long_abun, taxa(object, "tbl"), by = "taxa_id")
    }

    # LOOP OVER THE EXPRESSIONS AND GROUPING VARIABLES
    #----------------------------------------------------------------------------#
    for (i in seq_along(expressions)) {
      for(var in vars){
        
        expr_vars <- all.vars(expressions[[i]])
        expr_name <- names(expressions)[i]
        new_columns_name <- paste("color", var, expr_name, sep = "_")

        if (any(expr_vars %in% c("abun", "rela", "norm"))) {
          
          var_sorted <- long_abun %>%
            group_by(!!!rlang::syms(var)) %>%
            reframe(!!expressions[[i]]) %>%
            as.data.frame()
          
          var_sorted <- var_sorted[order(var_sorted[, ncol(var_sorted)], decreasing = decreasing), ]
          var_sorted <- pull(var_sorted, {{var}})
          
        } else {
          
          var_sorted <- taxa_info_mutated %>%
            group_by(!!!rlang::syms(var)) %>%
            reframe(!!expressions[[i]]) %>%
            as.data.frame()
          
          var_sorted <- var_sorted[order(var_sorted[, ncol(var_sorted)], decreasing = decreasing), ]
          var_sorted <- pull(var_sorted, {{var}})
          
        }
        
        distinct_color_var <- min(distinct_colors, length(var_sorted))
        palette_colors_var_expri <- colormap_categories(var_sorted,
                                                        distinct_colors = distinct_color_var,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        taxa_info_mutated <- taxa_info_mutated %>%
          mutate(!!new_columns_name := palette_colors_var_expri[taxa_info_mutated[[{{var}}]]])
        
    }} # end i,var
    
    taxa(object) <- taxa_info_mutated %>%
      dplyr::select(-tidyselect::any_of("comm_id")) %>%
      column_to_rownames("taxa_id")
    return(object)
  })


setMethod("set_categorical_color", signature = "mgnetList", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                                                     distinct_colors = 20, colorspace = "pretty",
                                                                     alpha = 1, extraColor = "#FFFFFF"){
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  expr_names <- names(expressions)
  
  # Check if any of the expressions are unnamed
  if (any(nzchar(expr_names) == 0)) {
    stop("All the expressions must be named.")
  }
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Check all expressions are named
  if(any(nzchar(expr_names)==0) | any(duplicated(expr_names))){
    stop("All the expressions must have an unique name")
  }
  
  # Capture required keys from expressions
  keys_required <- expressions %>%
    purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
    unlist() %>%
    unique()
  
  # Store needed abundances keys
  needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
  
  # Check the variables needed
  sapply(object, \(x) validate_required_variables(x, keys_required, "taxa"))
  
  if (!is.character(vars) || "sample_id" %in% vars) {
    stop("Error: 'vars' must be a character vector and cannot include 'sample_id'.")
  }
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)

  # Check the variables needed
  lapply(object, \(x){
    validate_required_variables(x, keys_required, "taxa")})
  
  # Generate new column names by combining variable elements with expression names
  combinations <- expand.grid(var = vars, expr = expr_names)
  new_column_names <- apply(combinations, 1, function(x) paste("color", x["var"], x["expr"], sep = "_"))

  # Check if these new column names already exist in the taxa
  existing_columns <- colnames(taxa(object, .fmt = "tbl"))
  if (any(new_column_names %in% existing_columns)) {
    conflicting_names <- new_column_names[new_column_names %in% existing_columns]
    warning("The following column names in taxa data will be overwritten: ", paste(conflicting_names, collapse = ", "))
  }
  
  if(distinct_colors > 99 || distinct_colors < 0){
    stop("Error: Ddstinct_colors must be in range [0,99].")
  }
  
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Error: alpha must be in the range [0,1].")
  }
  
  if(!is.character(colorspace) && !is.list(colorspace)) {
    stop("Error: colorspace must be a character or a list.")
  }
  
  if(is.character(colorspace)){
    if(!colorspace %in% c("pretty", "pretty_dark", "rainbow", "pastels")){
      stop("Error: when colorspace is a character, it must be one of 'pretty', 'pretty_dark', 'rainbow', 'pastels'.")
    }
  }
  
  if(is.list(colorspace) && !all(names(colorspace) %in% c("h", "s", "l"))) {
    stop("Error: when colorspace is a list, it must contain 'h', 's', and 'l' elements as required by qualpalr.")
  }
  
  if(!is.logical(decreasing)) stop("Error: decreasing must be logical.")
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # MERGED
  #----------------------------------------------------------------------------#
  #----------------------------------------------------------------------------#
  if(mode == "merged"){
    
    # Compile list of valid variables including "taxa_id" and additional taxa variables
    valid_vars <- c("taxa_id", unique(unlist(taxa_vars(object))))
    
    # Check if all elements of 'var' are in 'valid_vars'
    if (!all(vars %in% valid_vars)) {
      missing_vars <- vars[!vars %in% valid_vars]
      stop("Error: The following grouping 'vars' elements are not valid or missing in the taxa data: ", 
           paste(missing_vars, collapse = ", "))
    }
    
    # CREATE THE BASE FOR THE SOLUTION
    #----------------------------------------------------------------------------#
    taxa_info_mutated_merged <- object %>%
      purrr::map(\(x) {
        if(length(taxa(x)) != 0){
          taxa(x, .fmt = "tbl")
        } else {
          tibble::tibble(taxa_id = taxa_id(x))
        }
      }) %>%
      purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
      purrr::list_rbind()
    
    # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
    #----------------------------------------------------------------------------#
    if(length(needed_abundance_keys)!=0){
      
      long_abun_merged <- purrr::map(object, \(mgnet_obj){
        
        long_abun <- tidyr::expand_grid(sample_id = sample_id(mgnet_obj),
                                        taxa_id = taxa_id(mgnet_obj))
        
        for(abundance_key in needed_abundance_keys){
          long_abun <- long_abun %>%
            dplyr::left_join(methods::slot(mgnet_obj, abundance_key) %>%
                               tibble::as_tibble(rownames = "sample_id") %>%
                               tidyr::pivot_longer(-sample_id, 
                                                   names_to = "taxa_id", 
                                                   values_to = abundance_key),
                             by = c("sample_id", "taxa_id"))
        }
        
        if(length(taxa(mgnet_obj))!=0) long_abun <- dplyr::left_join(long_abun, taxa(mgnet_obj, "tbl"), by = "taxa_id")
        return(long_abun)
      }) %>%
        purrr::imap(\(x, y) tibble::add_column(x, mgnet = y, .before = 1)) %>%
        purrr::list_rbind()
    }
    
    # LOOP OVER THE EXPRESSIONS AND GROUPING VARIABLES
    #----------------------------------------------------------------------------#
    for (i in seq_along(expressions)) {
      for(var in vars){
        
        expr_vars <- all.vars(expressions[[i]])
        expr_name <- names(expressions)[i]
        new_columns_name <- paste("color", var, expr_name, sep = "_")
        
        if (any(expr_vars %in% c("abun", "rela", "norm"))) {
          
          var_sorted <- long_abun_merged %>%
            group_by(!!!rlang::syms(var)) %>%
            reframe(!!expressions[[i]]) %>%
            as.data.frame()
          
          var_sorted <- var_sorted[order(var_sorted[, ncol(var_sorted)], decreasing = decreasing), ]
          var_sorted <- pull(var_sorted, {{var}})
          
        } else {
          
          var_sorted <- taxa_info_mutated_merged %>%
            group_by(!!!rlang::syms(var)) %>%
            reframe(!!expressions[[i]]) %>%
            as.data.frame()
          
          var_sorted <- var_sorted[order(var_sorted[, ncol(var_sorted)], decreasing = decreasing), ]
          var_sorted <- pull(var_sorted, {{var}})
          
        }
        
        distinct_color_var <- min(distinct_colors, length(var_sorted))
        palette_colors_var_expri <- colormap_categories(var_sorted,
                                                        distinct_colors = distinct_color_var,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        taxa_info_mutated_merged <- taxa_info_mutated_merged %>%
          mutate(!!new_columns_name := palette_colors_var_expri[taxa_info_mutated_merged[[{{var}}]]])
        
      }} # end i,var
    
    taxa_info_mutated_splitted <- taxa_info_mutated_merged %>%
      base::split(.[, "mgnet"]) %>%
      purrr::imap(\(x,y){
        dplyr::arrange(x, match(taxa_id, taxa_id(object[[y]])))
      }) %>%
      purrr::map(\(x){
        x %>% dplyr::select(-any_of(c("mgnet", "comm_id"))) %>%
          tibble::column_to_rownames("taxa_id")
      })
    
    taxa_info_mutated_splitted <- taxa_info_mutated_splitted[names(object)]
    taxa(object) <- taxa_info_mutated_splitted
    return(object)
    
  } else if(mode == "separate"){
    
    # Check that 'by' is present in all mgnets
    for(x in object@mgnets) {
      if(!all((vars %in% c("taxa_id", taxa_vars(x))))) {
        stop(paste("Error: 'by' variable", vars, "is not present in all mgnet objects' taxa metadata."))
      }
    }
      
    # LOOP OVER THE EXPRESSIONS AND GROUPING VARIABLES
    #----------------------------------------------------------------------------#
    for (i in seq_along(expressions)) {
      for(var in vars){
        
        expr_vars <- all.vars(expressions[[i]])
        expr_name <- names(expressions)[i]
        new_columns_name <- paste("color", var, expr_name, sep = "_")
        
        var_sorted_separated <- top_taxa(object, expression = !!expressions[[i]], by = var, 
                                         decreasing = decreasing, mode = "separate")
        
        var_sorted_union_distinct <- sapply(var_sorted_separated, \(x) x[1:distinct_colors],
                                            simplify = FALSE, USE.NAMES = TRUE) %>%
          unlist() %>% unique()
        
        distinct_color_union <- length(var_sorted_union_distinct)
        var_sorted_all <- unique(c(var_sorted_union_distinct, unlist(var_sorted_separated)))
        
        palette_colors_var_expri <- colormap_categories(var_sorted_all,
                                                        distinct_colors = distinct_color_union,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        taxa(object) <- sapply(object, \(x){
          x <- taxa(x, "tbl")
          x <- x %>% 
            mutate(!!new_columns_name := palette_colors_var_expri[x[[{{var}}]]]) %>%
            dplyr::select(-any_of("comm_id")) %>%
            tibble::column_to_rownames("taxa_id")
          return(x)
          
        }, simplify = FALSE, USE.NAMES = TRUE)

    }} # end i var

    return(object)
  } # end separate
  
  
})