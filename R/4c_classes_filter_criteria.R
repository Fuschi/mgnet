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
#' @param log_abundance_criteria A list of functions to be applied to the log-transformed abundance data
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
           function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL, condition = "AND") {
             standardGeneric("filter_criteria_sample")
           })

setMethod("filter_criteria_sample", "mgnet",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL, condition = "AND") {
            
            # Checks
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(log_abundance_criteria)) return(object)
            
            # Check if all criteria are list, function, or NULL
            check_criteria <- function(crit) {
              if (is.null(crit)) return(TRUE)
              if (is.function(crit)) return(TRUE)
              if (is.list(crit) && all(sapply(crit, is.function))) return(TRUE)
              FALSE
            }
            if (!all(sapply(list(abundance_criteria, relative_criteria, log_abundance_criteria), check_criteria))) {
              stop("All criteria must be either a list of functions, a single function, or NULL.")
            }
            
            # Convert single function to list
            to_list <- function(x) if(is.function(x)) list(x) else x
            abundance_criteria <- to_list(abundance_criteria)
            relative_criteria <- to_list(relative_criteria)
            log_abundance_criteria <- to_list(log_abundance_criteria)
            
            apply_criteria <- function(data, criteria, condition){
              
              if(length(data)==0) return(data)
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
            # relative_pass <- if(length(relative(object))!=0){
            #   apply_criteria(relative(object), relative_criteria, condition)
            # } else {
            #   rep(TRUE, nsample(object))
            # }
            log_abundance_pass <- if(length(log_abundance(object))!=0){
              apply_criteria(log_abundance(object), log_abundance_criteria, condition)
            } else {
              rep(TRUE, nsample(object))
            }
            # Combine results based on condition
            final_pass <- if (condition == "AND") {
              abundance_pass & relative_pass & log_abundance_pass
            } else {
              abundance_pass | relative_pass | log_abundance_pass
            }
            
            subsetted_object <- object[final_pass, ]
            return(subsetted_object)
          })

setMethod("filter_criteria_sample", "mgnetList",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL, condition = "AND") {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              filter_criteria_sample(mgnet_obj, abundance_criteria, relative_criteria, log_abundance_criteria, condition)
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
#' @param log_abundance_criteria A list of functions to be applied to the log-transformed abundance data matrix column-wise.
#'        Each function takes a single column of log-transformed abundance data as input and returns a logical value.
#' @param condition A character string specifying the logical condition to apply across
#'        the criteria. Options are "AND" (default) and "OR". "AND" requires a taxa to meet all specified
#'        criteria to be included, while "OR" requires a taxa to meet at least one criterion.
#' @param trim A parameter to handle taxa not meeting criteria: "yes" to remove them,
#'        "no" to set their abundance to zero, or "aggregate" to combine them into a 'merged' category.
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
#' @seealso \link{filter_criteria_sample} for sample-based criteria application.
#' @export
#' @name filter_criteria_taxa
#' @aliases filter_criteria_taxa,mgnet-method filter_criteria_taxa,mgnetList-method
setGeneric("filter_criteria_taxa",
           function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL,
                    condition = "AND", trim = "yes", exclude_absent_taxa = TRUE) {
             standardGeneric("filter_criteria_taxa")
           })

setMethod("filter_criteria_taxa", "mgnet",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL,
                   condition = "AND", trim = "yes", exclude_absent_taxa = TRUE) {
            
            # Checks
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(log_abundance_criteria)) return(object)
            
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
            if (!all(sapply(list(abundance_criteria, relative_criteria, log_abundance_criteria), check_criteria))) {
              stop("All criteria must be either a list of functions, a single function, or NULL.")
            }
            
            # Convert single function to list
            to_list <- function(x) if(is.function(x)) list(x) else x
            abundance_criteria <- to_list(abundance_criteria)
            relative_criteria <- to_list(relative_criteria)
            log_abundance_criteria <- to_list(log_abundance_criteria)
            
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
            relative_pass <- if(length(relative(object))!=0){
              apply_criteria(relative(object), relative_criteria, condition)
            } else {
              rep(TRUE, ntaxa(object))
            }
            log_abundance_pass <- if(length(log_abundance(object))!=0){
              apply_criteria(log_abundance(object), log_abundance_criteria, condition)
            } else {
              rep(TRUE, ntaxa(object))
            }
            # Combine results based on condition
            final_pass <- if (condition == "AND") {
              abundance_pass & relative_pass & log_abundance_pass
            } else {
              abundance_pass | relative_pass | log_abundance_pass
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
              
              if(length(object@log_abundance)!=0){
                log_abundance.new <- object@log_abundance
                sample_min <- apply(object@log_abundance, 1, min)
                log_abundance.new[,!final_pass] <- sample_min # This could be an error (see later)
              } else {
                log_abundance.new <- object@log_abundance
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
                           log_abundance=log_abundance.new,
                           network=network.new,
                           community=community.new))
              # aggregate
              #-----------------------------------------#  
            } else if (trim=="aggregate") {
              
              # data
              if(length(object@abundance)!=0){
                abundance.new<-object@abundance[,final_pass,drop=F]
                if("merged"%in%colnames(abundance.new)){
                  abundance.new[,"merged"] <- abundance.new[,"merged"] + rowSums(object@abundance[,!final_pass,drop=F])
                } else {
                  abundance.new <- cbind(abundance.new, "merged"=rowSums(object@abundance[,!final_pass,drop=F]))
                }
              } else {
                abundance.new<-object@abundance
              }
              # taxa
              if(length(object@lineage)!=0){
                lineage.new<-object@lineage[final_pass,,drop=F]
                if(!("merged"%in%rownames(lineage.new))){
                  lineage.new <- rbind(lineage.new,"merged"=rep("merged",length(ranks(object))))
                }
              } else {
                lineage.new<-object@lineage
              }
              # info_taxa
              if(length(object@info_taxa)!=0){
                info_lineage.new<-object@info_taxa[final_pass,,drop=F]
                if(!("merged"%in%rownames(info_lineage.new))){
                  info_lineage.new <- rbind(info_lineage.new,"merged"=rep("merged",ncol(object@info_taxa)))
                }
              } else {
                info_lineage.new<-object@info_taxa
              }
              # log_data
              if(length(object@log_abundance)!=0){
                log_abundance.new<-object@log_abundance[,final_pass,drop=F]
                if("merged"%in%colnames(log_abundance.new)){
                  log_abundance.new[,"merged"] <- log_abundance.new[,"merged"] + rowSums(object@log_abundance[,!final_pass,drop=F])
                } else {
                  log_abundance.new <- cbind(log_abundance.new, "merged"=rowSums(object@log_abundance[,!final_pass,drop=F]))
                }
              } else {
                log_abundance.new<-object@log_abundance
              }
              # netw
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,final_pass)
                if(!("merged"%in%taxa_id(object))){
                  network.new <- igraph::add_vertices(network.new, nv=1, attr=list("name"="merged"))
                }
              } else {
                network.new<-object@network
              }
              # comm
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(final_pass)) final_pass <- which(taxa_id(object)%in%final_pass)
                community.new$membership <- object@community$membership[final_pass]
                if(!("merged"%in%taxa_id(object))){
                  community.new$membership <- c(community.new$membership,"merged"=0)
                }
                community.new$vcount <- length(community.new$membership)
                community.new$modularity <- NA
              } else {
                community.new <- object@community
              }
              
              return(mgnet(abundance=abundance.new,
                           info_sample=object@info_sample,
                           lineage=lineage.new,
                           info_taxa=info_lineage.new,
                           log_abundance=log_abundance.new,
                           network=network.new,
                           community=community.new))
            }
            
            subsetted_object <- object[ , final_pass, drop = FALSE]
            return(subsetted_object)
          })

setMethod("filter_criteria_taxa", "mgnetList",
          function(object, abundance_criteria = NULL, relative_criteria = NULL, log_abundance_criteria = NULL, 
                   condition = "AND", trim = "yes", exclude_absent_taxa = TRUE ) {
            # Apply the filtering criteria to each mgnet object in the list
            object@mgnets <- lapply(object@mgnets, function(mgnet_obj) {
              filter_criteria_taxa(mgnet_obj, abundance_criteria, relative_criteria, log_abundance_criteria, 
                                   condition, trim, exclude_absent_taxa)
            })
            return(object)
          })