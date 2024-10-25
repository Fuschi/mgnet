# SET CATEGORICAL COLOR
#------------------------------------------------------------------------------#
#' Set Categorical Colors in `mgnet` or `mgnetList` Objects
#' 
#' Documentation Not Completed
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
#' @param n Integer; the maximum number of distinct colors to generate for the categories. Defaults to 20.
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
#' @importFrom dplyr mutate group_by arrange
#' @importFrom tidyr expand_grid
#' @importFrom purrr map_chr
#' @importFrom methods slot
#' @importFrom tibble as_tibble
#' @export
#' @importFrom dplyr mutate %>%
#' @aliases set_categorical_color,mgnet-method set_categorical_color,mgnetList-method
setGeneric("set_categorical_color", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                             n = 20, colorspace = "pretty",
                                             alpha = 1, extraColor = "#FFFFFF"
                                             ) standardGeneric("set_categorical_color"))

setMethod("set_categorical_color", signature = "mgnet", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                                                 n = 20, colorspace = "pretty",
                                                                 alpha = 1, extraColor = "#FFFFFF"){

    # CHECKS
    #----------------------------------------------------------------------------#
  
    # Capture all the expressions provided
    expressions <- rlang::enquos(...)
    expr_names <- names(expressions)
    
    # Check all expressions have unique names
    if(any(nzchar(expr_names)==0) | any(duplicated(expr_names)))stop("All the expressions must have an unique name")

    # Check the reserved keywords
    check_reserved_keywords(expressions)

    # Check the variables needed
    needed_keys <- validate_required_variables(object, expressions, "taxa", FALSE)
    vars <- validate_required_groups(object, vars, "taxa")

    # Forbidden functions and disallowed variables
    check_forbidden_expressions(expressions)
    
    # Generate new column names by combining variable elements with expression names
    new_column_names <- purrr::map2_chr(vars, expr_names, ~ paste("color", .x, .y, sep = "_"))
    
    # Check if these new column names already exist in the taxa
    existing_columns <- colnames(taxa(object, .fmt = "tbl"))
    if (any(new_column_names %in% existing_columns)) {
      conflicting_names <- new_column_names[new_column_names %in% existing_columns]
      warning("The following column names in taxa data will be overwritten: ", paste(conflicting_names, collapse = ", "))
    }
    
    if(!is.logical(decreasing)) stop("Error: decreasing must be logical.")
    
    # END CHECKS
    #----------------------------------------------------------------------------#

    # CREATE THE BASE FOR THE SOLUTION
    #----------------------------------------------------------------------------#
    taxa_mutated <- initialize_taxa(object)
    long_abun <- long_abundance_join(object, needed_keys$abundance) %>%
      left_join(taxa_mutated, by = "taxa_id")

    # LOOP OVER THE EXPRESSIONS AND GROUPING VARIABLES
    #----------------------------------------------------------------------------#
    for (i in seq_along(expressions)) {
      for(var in vars){

        result <- apply_expr_var_sort(taxa_mutated, long_abun, expressions[[i]], var, decreasing, n)
        
        palette_colors_var_expri <- colormap_categories(as.character(result$var_sorted),
                                                        distinct_colors = result$n,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        new_columns_name <- paste("color", var, names(expressions)[i], sep = "_")
        taxa_mutated <- taxa_mutated %>%
          mutate(!!new_columns_name := palette_colors_var_expri[taxa_mutated[[{{var}}]]])
        
    }} # end i,var
    
    taxa(object) <- taxa_mutated %>%
      dplyr::select(-tidyselect::any_of("comm_id")) %>%
      column_to_rownames("taxa_id")
    
    return(object)
  })


setMethod("set_categorical_color", signature = "mgnetList", function(object, ..., vars, decreasing = TRUE, mode = "merged",
                                                                     n = 20, colorspace = "pretty",
                                                                     alpha = 1, extraColor = "#FFFFFF"){
  
  # CHECKS
  #----------------------------------------------------------------------------#
  
  # Capture all the expressions provided
  expressions <- rlang::enquos(...)
  expr_names <- names(expressions)
  
  # Check all expressions have unique names
  if(any(nzchar(expr_names)==0) | any(duplicated(expr_names)))stop("All the expressions must have an unique name")
  
  # Check the reserved keywords
  check_reserved_keywords(expressions)
  
  # Check the variables needed
  needed_keys <- validate_required_variables(object, expressions, "taxa", FALSE)
  vars <- validate_required_groups(object, vars, "taxa")
  
  # Forbidden functions and disallowed variables
  check_forbidden_expressions(expressions)
  
  # Generate new column names by combining variable elements with expression names
  new_column_names <- purrr::map2_chr(vars, expr_names, ~ paste("color", .x, .y, sep = "_"))
  
  # Check if these new column names already exist in the taxa
  existing_columns <- colnames(taxa(object, .fmt = "tbl"))
  if (any(new_column_names %in% existing_columns)) {
    conflicting_names <- new_column_names[new_column_names %in% existing_columns]
    warning("The following column names in taxa data will be overwritten: ", paste(conflicting_names, collapse = ", "))
  }
  
  if(!is.logical(decreasing)) stop("Error: decreasing must be logical.")
  
  # END CHECKS
  #----------------------------------------------------------------------------#
  
  # CREATE THE BASE FOR THE SOLUTION
  #----------------------------------------------------------------------------#
  taxa_mutated_merged <- initialize_taxa(object)
  long_abun_merged <- long_abundance_join(object, needed_keys$abundance) %>%
    left_join(taxa_mutated_merged, by = c("mgnet", "taxa_id"))
  
  # MERGED
  #----------------------------------------------------------------------------#
  if(mode == "merged"){
    
    for (i in seq_along(expressions)) {
      for(var in vars){
        
        result <- apply_expr_var_sort(taxa_mutated_merged, long_abun_merged, expressions[[i]], var, n, decreasing)

        palette_colors_var_expri <- colormap_categories(as.character(result$var_sorted[[var]]),
                                                        distinct_colors = result$n,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        new_columns_name <- paste("color", var, names(expressions)[i], sep = "_")
        taxa_mutated_merged <- taxa_mutated_merged %>%
          mutate(!!new_columns_name := palette_colors_var_expri[taxa_mutated_merged[[{{var}}]]])
        
      }} # end i,var
    
    taxa(object) <- split_arrange_merged_taxa(taxa_mutated_merged, object)
    return(object)
    
    
  # SEPARATE
  #----------------------------------------------------------------------------#
  } else if(mode == "separate"){
    
    # Check that 'by' is present in all mgnets
    for(x in object@mgnets) {
      if(!all((vars %in% taxa_vars(x)))) {
        stop(paste("Error: 'by' variable", vars, "is not present in all mgnet objects' taxa metadata."))
      }
    }
    
    for (i in seq_along(expressions)) {
      for(var in vars){
        
        result <- apply_expr_var_sort(taxa_mutated_merged, long_abun_merged, expressions[[i]], c("mgnet", var), n, decreasing)
        
        if(is.numeric(result$var_sorted[["_internal_"]])){
          
          categories <- result$var_sorted %>%
            group_by(mgnet) %>%
            arrange(if(decreasing) desc(`_internal_`) else `_internal_`) %>%
            slice_head(n = result$n) %>%
            ungroup() %>%
            pull(var)
          
          n = length(categories)
          
        } else {
          
          categories <- result$var_sorted %>%
            filter(`_internal_`) %>%
            pull(var) %>%
            unique()
          n <- length(categories)
          
        }
        
        categories <- unique(c(categories, result$var_sorted[[var]]))
        palette_colors_var_expri <- colormap_categories(as.character(categories),
                                                        distinct_colors = n,
                                                        colorspace = colorspace,
                                                        extraColor = extraColor,
                                                        alpha = alpha)
        
        new_columns_name <- paste("color", var, names(expressions)[i], sep = "_")
        taxa_mutated_merged <- taxa_mutated_merged %>%
          mutate(!!new_columns_name := palette_colors_var_expri[taxa_mutated_merged[[{{var}}]]])
        
      }} # end i,var
    
    taxa(object) <- split_arrange_merged_taxa(taxa_mutated_merged, object)
    return(object)
      
  
  } # end separate
  
  
})


# SET COMMUNITIES COLOR
#------------------------------------------------------------------------------#
#' Set Community Colors in `mgnet` Objects
#'
#' This function applies a color scheme to taxa based on community membership, which is determined by 
#' the sizes of communities identified in a metagenomic network. Each community can be assigned a unique 
#' color, while smaller communities and isolated nodes are assigned default colors to maintain visual 
#' clarity.
#'
#' @description
#' `set_community_color` assigns colors to taxa based on community memberships that are specified in the `comm_id` 
#' of the `mgnet` object. The function employs a color palette that is sensitive to community size, differentiating 
#' between larger and smaller communities.
#'
#' @param object An `mgnet` object that includes community data.
#' @param size An integer that specifies the minimum size a community must have to receive a unique color.
#'        Communities smaller than this size will be colored with `smaller_color` (default 1).
#' @param color_to The name of the column in the taxa data where the color values will be stored (default 'color_comm').
#' @param smaller_color A hexadecimal color code used for communities that do not exceed the `size` threshold,
#'        defaulting to a light gray ("#D3D3D3").
#' @param isolated_color A hexadecimal color code used for isolated nodes, often those not belonging to any community,
#'        defaulting to white ("#FFFFFF").
#' @param colorspace The qualitative color palette from which to generate the community colors. Options include
#'        'pretty', 'pretty_dark', 'rainbow', 'pastels', and others as supported by the `qualpalr` package. This
#'        parameter allows customization of the color scheme to fit different visualization needs.
#' @param alpha A numeric value between 0 and 1 indicating the opacity of the colors, where 1 is fully opaque and
#'        0 is completely transparent. This is useful for creating layered visual effects in plots.
#'
#' @details
#' The function modifies the provided `mgnet` object by adding or updating a column with community colors.
#' It relies on the `colormap_community` function to generate a suitable color palette based on community
#' sizes and the specified colorspace. The visualization can then be tailored to highlight community structures
#' within metagenomic networks effectively.
#'
#' @return The modified `mgnet` object with the added or updated community color data.
#'
#' @export
#' @seealso \link[colormap_community]{Generate Color Palette for Community Visualization}
#' @aliases set_community_color,mgnet-method set_community_color,mgnetList-method
setGeneric("set_community_color", function(object, size, color_to = "color_comm",
                                           smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                           colorspace = "pretty", alpha = 1) standardGeneric("set_community_color"))

setMethod("set_community_color", signature = "mgnet", function(object, size, color_to = "color_comm",
                                                               smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                                               colorspace = "pretty", alpha = 1){
  
  palette <- colormap_community(comm(object), sizes_threshold = size, 
                                smaller_color = smaller_color, isolated_color = isolated_color,
                                colorspace = colorspace, alpha = alpha)
  
  
  object <- mutate_taxa(object, !!color_to := palette[as.character(comm_id)])
  return(object)
})


setMethod("set_community_color", signature = "mgnetList", function(object, size, color_to = "color_comm",
                                                               smaller_color = "#D3D3D3", isolated_color = "#FFFFFF",
                                                               colorspace = "pretty", alpha = 1){
  
  palette <- colormap_communities(comm(object), sizes_threshold = size, 
                                  smaller_color = smaller_color, isolated_color = isolated_color,
                                  colorspace = colorspace, alpha = alpha)
  
  for(i in names(mgl)){
    object[[i]] <- mutate_taxa(object[[i]], !!color_to := palette[[i]][as.character(comm_id)])
  }
  
  return(object)
})

