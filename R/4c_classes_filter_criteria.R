# FILTER CRITERIA SAMPLE
#------------------------------------------------------------------------------#
#' Apply Criteria to Filter Samples Across an mgnetList
#'
#' This function allows for filtering samples across all `mgnet` objects contained within an `mgnetList`
#' based on specified criteria applied to abundance, relative abundance, and log-transformed
#' abundance data. Custom functions define the criteria for filtering, allowing for complex and
#' flexible sample selection processes across multiple datasets.
#'
#' @param object An `mgnetList` object.
#' @param abundance_criteria A list of functions to be applied to the abundance data matrix.
#'        Each function should take a single row of abundance data as input and return a logical
#'        value indicating whether the sample meets the specified criteria. Default is `NULL`.
#' @param relative_criteria A list of functions to be applied to the relative abundance data matrix.
#'        Each function should take a single row of relative abundance data as input and return
#'        a logical value. Default is `NULL`.
#' @param normalized_criteria A list of functions to be applied to the log-transformed abundance data
#'        matrix. Each function should take a single row of log-transformed abundance data as input
#'        and return a logical value. Default is `NULL`.
#' @param condition A character string specifying the logical condition to apply across
#'        the criteria. Options are "AND" (default) and "OR". "AND" requires a sample
#'        to meet all specified criteria to be included, while "OR" requires a sample
#'        to meet at least one criterion.
#'
#' @details
#' The `filter_criteria_sample` method for `mgnetList` iterates over each `mgnet` object in the list,
#' applying the specified filtering criteria to each. This method is useful for batch processing of
#' multiple datasets where the same filtering logic is desired across all datasets.
#'
#' @return A modified `mgnetList` object containing only the samples that meet all
#'         specified criteria across all provided data types. Each `mgnet` object within the list
#'         is filtered accordingly.
#'
#' @export
#' @seealso \link{filter_criteria_taxa} for taxa-based criteria application.
#' @name filter_criteria_sample
#' @aliases filter_criteria_sample,mgnet-method filter_criteria_sample,mgnetList-method
setGeneric("filter_criteria_sample",
           function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL, condition = "AND") {
             standardGeneric("filter_criteria_sample")
           })

setMethod("filter_criteria_sample", "mgnet",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL, condition = "AND") {
            
            # Checks
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(normalized_criteria)) return(object)
            
            # Check if all criteria are list, function, or NULL
            check_criteria <- function(crit) {
              if (is.null(crit)) return(TRUE)
              if (is.function(crit)) return(TRUE)
              if (is.list(crit) && all(sapply(crit, is.function))) return(TRUE)
              FALSE
            }
            if (!all(sapply(list(abundance_criteria, relative_criteria, normalized_criteria), check_criteria))) {
              stop("All criteria must be either a list of functions, a single function, or NULL.")
            }
            
            # Convert single function to list
            to_list <- function(x) if(is.function(x)) list(x) else x
            abundance_criteria <- to_list(abundance_criteria)
            relative_criteria <- to_list(relative_criteria)
            normalized_criteria <- to_list(normalized_criteria)
            
            apply_criteria <- function(data, criteria, condition){
              
              if(length(data)==0) stop(paste(data, "slot missing."))
              if (is.null(criteria)) return(rep(TRUE,nrow(data)))
              
              result_mat <- matrix(FALSE, nrow=nrow(data), ncol=length(criteria))
              for (r in 1:nrow(data)) {
                for (c in 1:length(criteria)) {
                  
                  result <- criteria[[c]](data[r,])
                  if (!is.logical(result) || length(result) != 1) {stop("Each criterion function must return a single logical value.")}
                  result_mat[r,c] <- result
                  
                }}
              
              if(condition=="AND"){
                pass <- as.logical(apply(result_mat,1,prod))
              } else {
                pass <- apply(result_mat,2,sum)>0
              }
              
              return(pass)
            }
            
            
            # Apply criteria functions to each data type
            abundance_pass <- if(length(abundance(object))!=0){
              apply_criteria(abundance(object), abundance_criteria, condition)
            } else {
              rep(TRUE, nsample(object))
            }
            relative_pass <- if(length(rel_abundance(object))!=0){
              apply_criteria(rel_abundance(object), relative_criteria, condition)
            } else {
              rep(TRUE, nsample(object))
            }
            norm_abundance_pass <- if(length(norm_abundance(object))!=0){
              apply_criteria(norm_abundance(object), normalized_criteria, condition)
            } else {
              rep(TRUE, nsample(object))
            }
            # Combine results based on condition
            final_pass <- if (condition == "AND") {
              abundance_pass & relative_pass & norm_abundance_pass
            } else {
              abundance_pass | relative_pass | norm_abundance_pass
            }
            
            subsetted_object <- object[final_pass, ]
            return(subsetted_object)
          })

setMethod("filter_criteria_sample", "mgnetList",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL, condition = "AND") {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              filter_criteria_sample(mgnet_obj, abundance_criteria, relative_criteria, normalized_criteria, condition)
            })
            return(object)
          })


# FILTER CRITERIA TAXA
#------------------------------------------------------------------------------#
#' Apply Criteria to Filter Taxa in mgnet Objects
#'
#' This function allows for filtering taxa in an `mgnet` object or across
#' an `mgnetList` based on specified criteria applied to abundance, relative abundance,
#' and log-transformed abundance data column-wise. Custom functions define the criteria for filtering,
#' allowing for complex and flexible taxa selection processes.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param abundance_criteria A list of functions to be applied to the abundance data matrix column-wise.
#'        Each function takes a single column of abundance data as input and returns a logical value
#'        indicating whether the taxa meets the filtering criteria.
#' @param relative_criteria A list of functions to be applied to the relative abundance data matrix column-wise.
#'        Each function takes a single column of relative abundance data as input and returns a logical value.
#' @param normalized_criteria A list of functions to be applied to the log-transformed abundance data matrix column-wise.
#'        Each function takes a single column of log-transformed abundance data as input and returns a logical value.
#' @param condition A character string specifying the logical condition to apply across
#'        the criteria. Options are "AND" (default) and "OR". "AND" requires a taxa to meet all specified
#'        criteria to be included, while "OR" requires a taxa to meet at least one criterion.
#' @param trim How to handle taxa not meeting criteria: "yes" to remove them,
#'        "no" to set their abundance to zero (min sample value for norm_abundance), 
#'        or "aggregate" to combine them into a single category.
#' @param aggregate_to Name for the aggregated taxa category when using `trim = "aggregate"`.
#'        Defaults to "aggregate".   
#' @param exclude_absent_taxa A logical value indicating whether to exclude taxa that are absent across all samples
#'        before applying the criteria. Defaults to TRUE.
#'
#' @details
#' The function evaluates each specified criterion against the data columns (taxa)
#' in the abundance, relative abundance, and log-transformed abundance matrices.
#' Taxa that meet all provided criteria across these datasets are retained in the filtered object.
#'
#' @return A filtered `mgnet` object containing only the taxa that meet all
#'         specified criteria, or an `mgnetList` with each constituent `mgnet`
#'         object filtered accordingly.
#'
#' @importFrom rlang !! :=
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#'
#' @seealso \link{filter_criteria_sample} for sample-based criteria application.
#' @export
#' @name filter_criteria_taxa
#' @aliases filter_criteria_taxa,mgnet-method filter_criteria_taxa,mgnetList-method
setGeneric("filter_criteria_taxa",
           function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL,
                    condition = "AND", trim = "yes", aggregate_to = NULL,
                    exclude_absent_taxa = TRUE) {
             standardGeneric("filter_criteria_taxa")
           })

setMethod("filter_criteria_taxa", "mgnet",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL,
                   condition = "AND", trim = "yes", aggregate_to = NULL,
                   exclude_absent_taxa = TRUE) {
            
            # Checks
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(normalized_criteria)) return(object)
            
            # Automatically remove taxa that are always equal to zero across all samples
            if(exclude_absent_taxa && length(abundance) != 0) {
              zero_taxa <- colSums(abundance(object)) != 0
              object <- object[, zero_taxa]
            }
            
            # Check if all criteria are list, function, or NULL
            check_criteria <- function(crit) {
              if (is.null(crit)) return(TRUE)
              if (is.function(crit)) return(TRUE)
              if (is.list(crit) && all(sapply(crit, is.function))) return(TRUE)
              FALSE
            }
            if (!all(sapply(list(abundance_criteria, relative_criteria, normalized_criteria), check_criteria))) {
              stop("All criteria must be either a list of functions, a single function, or NULL.")
            }
            
            if(trim == "aggregate" && is.null(aggregate_to)) stop("if you set trim as 'aggregate' the parameter aggregate_to must be set as character and indicate the new variable where the filtered taxa are aggregated")
            if(!is.null(aggregate_to)){
              if(!is.character(aggregate_to)) stop("aggregate_to must to be a character")
            }
            if(!is.null(aggregate_to) && trim != "aggregate"){warning("if the parameter trim is different from 'aggregate' the parameter aggregate to is ignored")}
            
            # Convert single function to list
            to_list <- function(x) if(is.function(x)) list(x) else x
            abundance_criteria <- to_list(abundance_criteria)
            relative_criteria <- to_list(relative_criteria)
            normalized_criteria <- to_list(normalized_criteria)
            
            apply_criteria <- function(data, criteria, condition){
              
              if(length(data)==0) return(data)
              
              if (is.null(criteria)) return(rep(TRUE,ncol(data)))
              
              result_mat <- matrix(FALSE, nrow=ncol(data), ncol=length(criteria))
              for (col in 1:ncol(data)) {
                for (cri in 1:length(criteria)) {
                  
                  result <- criteria[[cri]](data[,col])
                  if (is.na(result) || !is.logical(result) || length(result) != 1) {stop("Each criterion function must return a single logical value (TRUE | FALSE).")}
                  result_mat[col,cri] <- result
                  
                }}
              
              if(condition=="AND"){
                pass <- as.logical(apply(result_mat,1,prod))
              } else {
                pass <- apply(result_mat,1,sum)>0
              }
              return(pass)
            }
            
            
            # Apply criteria functions to each data type
            abundance_pass <- if(length(abundance(object))!=0){
              apply_criteria(abundance(object), abundance_criteria, condition)
            } else {
              rep(TRUE, ntaxa(object))
            }
            relative_pass <- if(length(rel_abundance(object))!=0){
              apply_criteria(rel_abundance(object), relative_criteria, condition)
            } else {
              rep(TRUE, ntaxa(object))
            }
            norm_abundance_pass <- if(length(norm_abundance(object))!=0){
              apply_criteria(norm_abundance(object), normalized_criteria, condition)
            } else {
              rep(TRUE, ntaxa(object))
            }
            # Combine results based on condition
            final_pass <- if (condition == "AND") {
              abundance_pass & relative_pass & norm_abundance_pass
            } else {
              abundance_pass | relative_pass | norm_abundance_pass
            }
            
            
            # yes
            #-----------------------------------------#
            if(trim=="yes"){
              
              return(object[ ,final_pass])
              
              # no
              #-----------------------------------------#
            } else if (trim=="no") {
              
              if(length(object@abundance)!=0){
                abundance.new <- object@abundance
                abundance.new[,!final_pass] <- 0
              } else {
                abundance.new <- object@abundance
              }
              
              if(length(object@norm_abundance)!=0){
                norm_abundance.new <- object@norm_abundance
                sample_min <- apply(object@norm_abundance, 1, min)
                norm_abundance.new[,!final_pass] <- sample_min # This could be an error (see later)
              } else {
                norm_abundance.new <- object@norm_abundance
              }
              
              if(length(object@network)!=0){
                sub.network <- subgraph(network(object), taxa_id(object)[final_pass])
                preserved.edges <- E(object@network)%in%E(sub.network)
                network.new <- subgraph.edges(graph=network(object), 
                                              eids=E(network(object))[preserved.edges],
                                              delete.vertices = FALSE)
              } else {
                network.new<-object@network
              }
              
              if(length(object@community)!=0){
                community.new <- object@community
                community.new$membership[!final_pass] <- 0
                community.new$modularity <- NA
              } else {
                community.new <- object@community
              }
              
              return(mgnet(abundance=abundance.new,
                           info_sample=object@info_sample,
                           lineage=object@lineage,
                           info_taxa=object@info_taxa,
                           norm_abundance=norm_abundance.new,
                           network=network.new,
                           community=community.new))
              # aggregate
              #-----------------------------------------#  
            } else if (trim=="aggregate") {
              
              # abundance
              if(length(object@abundance)!=0){
                abundance.new<-object@abundance[,final_pass,drop=F]
                if(aggregate_to%in%colnames(abundance.new)){
                  abundance.new[,aggregate_to] <- abundance.new[,aggregate_to] + rowSums(object@abundance[,!final_pass,drop=F])
                } else {
                  abundance.new <- as.data.frame(abundance.new) %>%
                    mutate(!!aggregate_to := rowSums(object@abundance[,!final_pass,drop=F]) ) %>%
                    as.matrix
                }
              } else {
                abundance.new<-object@abundance
              }
              # rel_abundance
              if(length(object@rel_abundance)!=0){
                rel_abundance.new<-object@rel_abundance[,final_pass,drop=F]
                if(aggregate_to%in%colnames(rel_abundance.new)){
                  rel_abundance.new[,aggregate_to] <- rel_abundance.new[,aggregate_to] + rowSums(object@rel_abundance[,!final_pass,drop=F])
                } else {
                  rel_abundance.new <- as.data.frame(rel_abundance.new) %>%
                    mutate(!!aggregate_to := rowSums(object@rel_abundance[,!final_pass,drop=F]) ) %>%
                    as.matrix
                }
              } else {
                rel_abundance.new<-object@rel_abundance
              }
              # norm_data
              if(length(object@norm_abundance)!=0){
                norm_abundance.new<-object@norm_abundance.new[,final_pass,drop=F]
                if(aggregate_to%in%colnames(norm_abundance.new)){
                  norm_abundance.new[,aggregate_to] <- norm_abundance.new[,aggregate_to] + rowSums(object@norm_abundance[,!final_pass,drop=F])
                } else {
                  norm_abundance.new <- as.data.frame(norm_abundance.new) %>%
                    mutate(!!aggregate_to := rowSums(object@norm_abundance[,!final_pass,drop=F]) ) %>%
                    as.matrix
                }
              } else {
                norm_abundance.new<-object@norm_abundance
              }
              # lineage
              if(length(object@lineage)!=0){
                lineage.new<-object@lineage[final_pass,,drop=F]
                if(!(aggregate_to%in%rownames(lineage.new))){
                  lineage.new <- rbind(lineage.new,
                                       matrix(aggregate_to, nrow = 1, ncol = length(ranks(object)),
                                              dimnames = list(aggregate_to, ranks(object)))) 
                }
              } else {
                lineage.new<-object@lineage
              }
              # info_taxa
              if(length(object@info_taxa)!=0){
                info_taxa.new<-object@info_taxa[final_pass,,drop=F]
                if(!(aggregate_to%in%rownames(info_taxa.new))){
                  info_taxa.new <- as.data.frame(info_taxa.new) %>%
                    mutate(!!aggregate_to := aggregate_to ) %>%
                    as.matrix
                }
              } else {
                info_taxa.new<-object@info_taxa
              }
              # network
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,final_pass)
                if(!(aggregate_to%in%taxa_id(object))){
                  network.new <- igraph::add_vertices(network.new, nv=1, attr=list("name"=aggregate_to))
                }
              } else {
                network.new<-object@network
              }
              # community
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(final_pass)) final_pass <- which(taxa_id(object)%in%final_pass)
                community.new$membership <- object@community$membership[final_pass]
                if(!(aggregate_to%in%taxa_id(object))){
                  community.new$membership <- c(community.new$membership,aggregate_to=0)
                }
                community.new$vcount <- length(community.new$membership)
                community.new$modularity <- NA
              } else {
                community.new <- object@community
              }
              
              return(mgnet(abundance = abundance.new, 
                           rel_abundance = rel_abundance.new,
                           norm_abundance = norm_abundance.new,
                           info_sample = object@info_sample,
                           lineage = lineage.new,
                           info_taxa = info_taxa.new,
                           network = network.new,
                           community = community.new))
            }

          })

setMethod("filter_criteria_taxa", "mgnetList",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, normalized_criteria = NULL, 
                   condition = "AND", trim = "yes", aggregate_to = NULL,
                   exclude_absent_taxa = TRUE ) {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              filter_criteria_taxa(mgnet_obj, abundance_criteria, relative_criteria, normalized_criteria, 
                                   condition, trim, aggregate_to, exclude_absent_taxa)
            })
            return(object)
          })