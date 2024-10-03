##' # SET CATEGORICAL COLOR#------------------------------------------------------------------------------#@export
##' @importFrom dplyr mutate %>%
##' @name set_categorical_color
##' @aliases set_categorical_color,mgnet-method set_categorical_color,mgnetList-method
# setGeneric("set_categorical_color", function(object, ..., by, distinct_colors, decreasing = FALSE) standardGeneric("set_categorical_color"))
# 
# setMethod("set_categorical_color", signature = "mgnet", function(object, ..., by, distinct_colors, decreasing = FALSE){
#     
#     # CHECKS
#     #----------------------------------------------------------------------------#
#     # Ensure there are taxa to process
#     if (ntaxa(object) == 0) {
#       stop("Error: No taxa available in the 'mgnet' object.")
#     }
#     
#     # Capture all the expressions provided
#     expressions <- rlang::enquos(...)
#     
#     # Check the reserved keywords
#     check_reserved_keywords(expressions)
#     
#     # Check all expressions are named
#     if(any(nzchar(names(expressions))==0)){
#       stop("All the expressions must be named")
#     }
#     
#     # Capture required keys from expressions
#     keys_required <- expressions %>%
#       purrr::map(~ all.vars(rlang::quo_get_expr(.x))) %>%
#       unlist() %>%
#       unique()
#     
#     # Store needed abundances keys
#     needed_abundance_keys <- intersect(keys_required, c("abun","rela","norm"))
#     
#     # Check the variables needed
#     validate_required_variables(object, keys_required, "taxa")
#     
#     # Validate the by argument
#     if (missing(by)) {
#       stop("Error: 'by' must be setted.")
#     }
#     
#     if (!is.character(by) || "sample_id" %in% by) {
#       stop("Error: 'by' must be a character vector and cannot include 'sample_id'.")
#     }
#     
#     # Forbidden functions and disallowed variables
#     check_forbidden_expressions(expressions)
#     
#     # END CHECKS
#     #----------------------------------------------------------------------------#
#     
#     # CREATE THE BASE FOR THE SOLUTION
#     #----------------------------------------------------------------------------#
#     # Initialize sample_info_mutated
#     if (length(taxa(object)) != 0) {
#       taxa_info_mutated <- taxa(object, .fmt = "tbl")
#     } else {
#       taxa_info_mutated <- tibble::tibble(taxa_id = taxa_id(object))
#     }
#     
#     # CALCULATE LONG ABUNDANCE DATA ONCE, if needed
#     #----------------------------------------------------------------------------#
#     if (length(needed_abundance_keys) > 0) {
#       # Create the grid of sample_id and taxa_id combinations
#       long_abun <- tidyr::expand_grid(
#         sample_id = sample_id(object),
#         taxa_id = taxa_id(object)
#       )
#       
#       # Left join all abundance variables needed
#       for (abundance_key in needed_abundance_keys) {
#         long_abun <- long_abun %>%
#           dplyr::left_join(
#             methods::slot(object, abundance_key) %>%
#               tibble::as_tibble(rownames = "sample_id") %>%
#               tidyr::pivot_longer(-sample_id,
#                                   names_to = "taxa_id",
#                                   values_to = abundance_key),
#             by = c("sample_id", "taxa_id")
#           )
#       }
#     }
#     
#     # LOOP OVER THE EXPRESSIONS
#     #----------------------------------------------------------------------------#
#     for (i in seq_along(expressions)) {
#       
#       expr_vars <- all.vars(expressions[[i]])
#       expr_name <- names(expressions)[i]
#       
#       if (any(expr_vars %in% c("abun", "rela", "norm"))) {
#         
#         # EXPRESSION WITH ABUNDANCES
#         #----------------------------------------------#
#         palette_tmp <- long_abun %>%
#           dplyr::left_join(taxa_info_mutated, by = "taxa_id") %>%
#           dplyr::group_by(!!!rlang::syms(by)) %>%
#           dplyr::summarise(expr_name := !!!rlang::eval_tidy(expressions[[i]])) %>%
#           dplyr::ungroup()# %>%
#           #dplyr::arrange(desc(!!rlang::sym(expr_name)))
#         return(palette_tmp)
#           #tidyr::unite(!!rlang::sym(expr_name), !!!rlang::syms(by), sep = "_", remove = FALSE) %>%
#           #dplyr::mutate(!!rlang::sym(expr_name) := colormap_categories(categories = !!rlang::sym(expr_name),
#           #                                                             distinct_colors = distinct_colors))
#         
#         return(palette_tmp)
#         
#       } else {
#         
#         # EXPRESSION WITHOUT ABBUNDANCES
#         #----------------------------------------------#
#         palette_tmp <- taxa_info_mutated %>%
#           dplyr::group_by(!!!rlang::syms(by)) %>%
#           dplyr::summarise(expr_name := !!!rlang::eval_tidy(expressions[i])) %>%
#           dplyr::ungroup() %>%
#           dplyr::arrange(desc(!!rlang::sym(expr_name))) %>%
#           mutate(expr_name := colormap_categories(!!rlang::sym(expr_name),
#                                                   distinct_colors = distinct_colors))
#         
#         taxa_info_mutated <- taxa_info_mutated %>%
#           left_join(palette_tmp, by = join_by(!!!rlang::syms(by)))
#           
#         
#         return(taxa_info_mutated)
#         
#       }
#     } # End of loop over expressions
#     
#   })