# REFINE SAMPLE 
#------------------------------------------------------------------------------#
#' Refines Sample Selection in `mgnet` or `mgnetList` Objects
#'
#' @description
#' This function refines the selection of samples within an `mgnet` object or across
#' `mgnet` objects within an `mgnetList`. It enables the specification of criteria
#' for which samples to retain, supporting logical vectors, numeric indices, sample names,
#' or custom functions for dynamic selection.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Criteria for sample selection: logical vectors, numeric indices, sample names,
#'        or functions returning such indices. These criteria determine which samples are retained.
#' @param condition Logical condition to combine multiple criteria: "AND" (default) or "OR".
#'        "AND" requires all criteria to be met, "OR" requires any criterion to be met.
#'
#' @return A modified `mgnet` object with only the selected samples, or an `mgnetList` with
#'         each contained `mgnet` object filtered accordingly.
#'
#' @details
#' This function allows for advanced sample selection by applying the specified criteria
#' directly to the sample data within an `mgnet` object. When provided with multiple criteria,
#' these are combined according to the specified logical condition. For `mgnetList` objects,
#' the function applies the selection process to each contained `mgnet` object independently,
#' facilitating bulk operations.
#'
#' @export
#' @aliases refine_sample,mgnet-method refine_sample,mgnetList-method
#' @seealso \link{mgnet}, \link{mgnetList}
setGeneric("refine_sample", function(object, ..., condition = "AND") standardGeneric("refine_sample"))

setMethod("refine_sample", "mgnet",
          function(object, ..., condition = "AND") {
            condition <- match.arg(condition, c("AND", "OR"))
            if (length(object@network)!=0) {
              warning("Sample subsetting removes network and community slots.")
            }
            
            condition <- match.arg(condition, c("AND","OR"))
            if(length(object@network)!=0) warning("sample subsetting is not defined for network, the resulting object will not have the slots netw and comm")
            
            IDX <- list(...)
            for(i in 1:length(IDX)){
              
              idx <- IDX[[i]]
              if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
              
              if(is.function(idx)){
                idx <- idx(object)
                if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
              }  
              
              
              if(is.numeric(idx)){
                #numeric
                if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>nsample(object)){
                  stop("numeric indices must be integers in range 1 to the maximum nuber of sample")
                }
                #transform to logical
                IDX[[i]] <- (1:nsample(object))%in%idx
                
              }else if(is.character(idx)){
                #character
                if( !any(idx%in%sample_id(object)) ){
                  stop("string indices must be a subset of sample_id of object")
                }
                #transform to logical
                IDX[[i]] <- sample_id(object)%in%idx
                
              }else if(is.logical(idx)){
                #logical
                if(length(idx)!=nsample(object)){
                  stop("logical indices must have the length equal to sample number")
                }
                IDX[[i]] <- idx
                
              } else {stop("the indices must character, numeric or logical vectors")}
            }
            
            
            IDX <- do.call(rbind,IDX)
            # Condition  
            if(condition=="AND"){
              IDX <- as.logical(apply(IDX,2,prod))
            } else {
              IDX <- apply(IDX,2,sum)>0
            }
            
            # Subset the mgnet object
            subsetted_object <- object[IDX, ]
            
            return(subsetted_object)
          }
)

setMethod("refine_sample", "mgnetList",
          function(object, ..., condition = "AND") {
            condition <- match.arg(condition, c("AND", "OR"))
            
            # Apply refine_taxa to each mgnet object in the mgnetList
            updated_mgnets <- lapply(object@mgnets, function(mgnetObj){
              refine_sample(object = mgnetObj, ... = ..., condition = condition)
            })
            
            # Update the mgnetList with the filtered mgnet objects
            object@mgnets <- updated_mgnets
            validObject(object)
            
            return(object)
          }
)


# REFINE TAXA 
#------------------------------------------------------------------------------#
#' Refines Taxa Selection in `mgnet` or `mgnetList` Objects
#'
#' @description
#' This function refines the selection of taxa within an `mgnet` object or across
#' `mgnet` objects within an `mgnetList`. It supports the specification of criteria
#' for which taxa to retain, using logical vectors, numeric indices, taxa names,
#' or functions that dynamically specify taxa selection criteria.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Criteria for taxa selection: logical vectors, numeric indices, taxa names,
#'        or functions returning such indices. These criteria determine which taxa are retained.
#' @param condition Logical condition to combine multiple criteria: "AND" (default) or "OR".
#'        "AND" requires all criteria to be met, "OR" requires any criterion to be met.
#' @param trim How to handle taxa not meeting criteria: "yes" to remove them,
#'        "no" to set their abundance to zero (min sample value for norm_abundance), 
#'        or "aggregate" to combine them into a single category.
#' @param aggregate_to Name for the aggregated taxa category when using `trim = "aggregate"`.
#'        Defaults to NULL.        
#'
#' @return A modified `mgnet` object with only the selected taxa, or an `mgnetList` with
#'         each contained `mgnet` object filtered accordingly.
#'
#' @details
#' This function provides a flexible approach to taxa refinement, allowing users to easily
#' include or exclude specific taxa based on custom logic. The `trim` parameter offers
#' additional control over how unselected taxa are treated, enabling tailored data cleaning
#' and preparation steps prior to further analysis.
#' 
#' @export
#' @aliases refine_taxa,mgnet-method  refine_taxa,mgnetList-method
#' @importFrom igraph subgraph add_vertices subgraph.edges E
#' @seealso \link{mgnet}, \link{mgnetList}
setGeneric("refine_taxa", function(object, ..., condition = "AND", 
                                   trim = "yes", aggregate_to = NULL) standardGeneric("refine_taxa"))

setMethod("refine_taxa", "mgnet",
          function(object, ..., condition = "AND",
                   trim = "yes", aggregate_to = NULL){
            
            condition <- match.arg(condition, c("AND","OR"))
            trim <- match.arg(trim, c("yes","no","aggregate"))
            
            if(trim == "aggregate" && is.null(aggregate_to)) stop("if you set trim as 'aggregate' the parameter aggregate_to must be set as character and indicate the new variable where the filtered taxa are aggregated")
            if(!is.null(aggregate_to)){
              if(!is.character(aggregate_to)) stop("aggregate_to must to be a character")
            }
            if(!is.null(aggregate_to) && trim != "aggregate"){warning("if the parameter trim is different from 'aggregate' the parameter aggregate to is ignored")}
            
            
            IDX <- list(...)
            for(i in 1:length(IDX)){
              
              idx <- IDX[[i]]
              if(!(is.vector(idx) | is.function(idx))) stop("all elements of ... must be vectors or functions")
              
              if(is.function(idx)){
                idx <- idx(object)
                if(!is.vector(idx))stop(paste("the function at position",i,"does not return a vector"))
              }  
              
              if(is.numeric(idx)){
                #numeric
                if(any(is.na(idx)) | any(idx<0) | any(round(idx)!=idx) | max(idx)>ntaxa(object)){
                  stop("numeric indices must be integers in range 1 to the maximum nuber of taxa")
                }
                #transform to logical
                IDX[[i]] <- (1:ntaxa(object))%in%idx
                
              }else if(is.character(idx)){
                #character
                if( !any(idx%in%taxa_id(object)) ){
                  stop("string indices must be a subset of taxa_id of object")
                }
                #transform to logical
                IDX[[i]] <- taxa_id(object)%in%idx
                
              }else if(is.logical(idx)){
                #logical
                if(length(idx)!=ntaxa(object)){
                  stop("logical indices must have the length equal to taxa number")
                }
                IDX[[i]] <- idx
                
              } else {stop("the indices must character, numeric or logical vectors")}
            }
            
            
            IDX <- do.call(rbind,IDX)
            # Condition  
            if(condition=="AND"){
              IDX <- as.logical(apply(IDX,2,prod))
            } else {
              IDX <- apply(IDX,2,sum)>0
            }
            
            # yes
            #-----------------------------------------#
            if(trim=="yes"){
              
              return(object[,IDX])
              
              # no
              #-----------------------------------------#
            } else if (trim=="no") {
              
              if(length(object@abundance)!=0){
                abundance.new <- object@abundance
                abundance.new[,!IDX] <- 0
              } else {
                abundance.new <- object@abundance
              }
              
              if(length(object@rel_abundance)!=0){
                rel_abundance.new <- object@rel_abundance
                rel_abundance.new[,!IDX] <- 0
              } else {
                rel_abundance.new <- object@rel_abundance
              }
              
              if(length(object@norm_abundance)!=0){
                norm_abundance.new <- object@norm_abundance
                sample_min <- apply(object@norm_abundance, 1, min)
                
                norm_abundance.new <- for(sample in seq_along(norm_abundance.new)){
                  norm_abundance.new[,!IDX] <- sample_min[sample]
                }
                
              } else {
                norm_abundance.new <- object@norm_abundance
              }
              
              if(length(object@network)!=0){
                sub.network <- subgraph(network(object), taxa_id(object)[IDX])
                preserved.edges <- E(object@network)%in%E(sub.network)
                network.new <- subgraph.edges(graph=network(object), 
                                              eids=E(network(object))[preserved.edges],
                                              delete.vertices = FALSE)
              } else {
                network.new<-object@network
              }
              
              if(length(object@community)!=0){
                community.new <- object@community
                community.new$membership[!IDX] <- 0
                community.new$modularity <- NA
              } else {
                community.new <- object@community
              }
              
              return(mgnet(abundance=abundance.new,
                           rel_abundance=rel_abundance.new,
                           norm_abundance=norm_abundance.new,
                           info_sample=object@info_sample,
                           lineage=object@lineage,
                           info_taxa=object@info_taxa,
                           network=network.new,
                           community=community.new))
              # aggregate
              #-----------------------------------------#  
            } else if (trim=="aggregate") {
              
              # abundance
              if(length(object@abundance)!=0){
                abundance.new<-object@abundance[,IDX,drop=F]
                if(aggregate_to%in%colnames(abundance.new)){
                  abundance.new[,aggregate_to] <- abundance.new[,aggregate_to] + rowSums(object@abundance[,!IDX,drop=F])
                } else {
                  abundance.new <- cbind(abundance.new, aggregate_to=rowSums(object@abundance[,!IDX,drop=F]))
                }
              } else {
                abundance.new<-object@abundance
              }
              # rel_abundance
              if(length(object@rel_abundance)!=0){
                rel_abundance.new<-object@rel_abundance[,IDX,drop=F]
                if(aggregate_to%in%colnames(rel_abundance.new)){
                  rel_abundance.new[,aggregate_to] <- rel_abundance.new[,aggregate_to] + rowSums(object@rel_abundance[,!IDX,drop=F])
                } else {
                  rel_abundance.new <- cbind(rel_abundance.new, aggregate_to=rowSums(object@rel_abundance[,!IDX,drop=F]))
                }
              } else {
                rel_abundance.new <- object@rel_abundance
              }
              # norm_abundance
              if(length(object@norm_abundance)!=0){
                norm_abundance.new<-object@norm_abundance[,IDX,drop=F]
                if(aggregate_to%in%colnames(norm_abundance.new)){
                  norm_abundance.new[,aggregate_to] <- norm_abundance.new[,aggregate_to] + rowSums(object@norm_abundance[,!IDX,drop=F])
                } else {
                  norm_abundance.new <- cbind(norm_abundance.new, aggregate_to=rowSums(object@norm_abundance[,!IDX,drop=F]))
                }
              } else {
                norm_abundance.new<-object@norm_abundance
              }
              # lineage
              if(length(object@lineage)!=0){
                lineage.new<-object@lineage[IDX,,drop=F]
                if(!(aggregate_to%in%rownames(lineage.new))){
                  lineage.new <- rbind(lineage.new,aggregate_to=rep(aggregate_to,length(ranks(object))))
                }
              } else {
                lineage.new<-object@lineage
              }
              # info_taxa
              if(length(object@info_taxa)!=0){
                info_lineage.new<-object@info_taxa[IDX,,drop=F]
                if(!(aggregate_to%in%rownames(info_lineage.new))){
                  info_lineage.new <- rbind(info_lineage.new,aggregate_to=rep(aggregate_to,ncol(object@info_taxa)))
                }
              } else {
                info_lineage.new<-object@info_taxa
              }
              # netw
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,IDX)
                if(!(aggregate_to%in%taxa_id(object))){
                  network.new <- igraph::add_vertices(network.new, nv=1, attr=list("name"=aggregate_to))
                }
              } else {
                network.new<-object@network
              }
              # comm
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(IDX)) IDX <- which(taxa_id(object)%in%IDX)
                community.new$membership <- object@community$membership[IDX]
                if(!(aggregate_to%in%taxa_id(object))){
                  community.new$membership <- c(community.new$membership,aggregate_to=0)
                }
                community.new$vcount <- length(community.new$membership)
                community.new$modularity <- NA
              } else {
                community.new <- object@community
              }
              
              return(mgnet(abundance=abundance.new,
                           rel_abundance=rel_abundance.new,
                           norm_abundance=norm_abundance.new,
                           info_sample=object@info_sample,
                           lineage=lineage.new,
                           info_taxa=info_lineage.new,
                           network=network.new,
                           community=community.new))
            }
            
          })


setMethod("refine_taxa", "mgnetList",
          function(object, ..., condition = "AND", 
                   trim = "yes", aggregate_to = NULL) {
            condition <- match.arg(condition, c("AND", "OR"))
            trim <- match.arg(trim, c("yes", "no", "aggregate"))
            
            # Apply refine_taxa to each mgnet object in the mgnetList
            updated_mgnets <- lapply(object@mgnets, function(mgnetObj){
              refine_taxa(object = mgnetObj, ... = ..., 
                          condition = condition, trim = trim)
            })
            
            # Update the mgnetList with the filtered mgnet objects
            object@mgnets <- updated_mgnets
            
            return(object)
          }
)