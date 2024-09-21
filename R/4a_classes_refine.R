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
            if (length(object@netw)!=0) {
              warning("Sample subsetting removes netw and comm slots.")
            }
            
            condition <- match.arg(condition, c("AND","OR"))
            if(length(object@netw)!=0) warning("sample subsetting is not defined for network, the resulting object will not have the slots netw and comm")
            
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
#'        "no" to set their abundance to zero (min sample value for norm), 
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
              
              if(length(object@abun)!=0){
                abun.new <- object@abun
                abun.new[,!IDX] <- 0
              } else {
                abun.new <- object@abun
              }
              
              if(length(object@rela)!=0){
                rela.new <- object@rela
                rela.new[,!IDX] <- 0
              } else {
                rela.new <- object@rela
              }
              
              if(length(object@norm)!=0){
                norm.new <- object@norm
                sample_min <- apply(object@norm, 1, min)
                
                norm.new <- for(sample in seq_along(norm.new)){
                  norm.new[,!IDX] <- sample_min[sample]
                }
                
              } else {
                norm.new <- object@norm
              }
              
              if(length(object@netw)!=0){
                sub.netw <- subgraph(netw(object), taxa_id(object)[IDX])
                preserved.edges <- E(object@netw)%in%E(sub.netw)
                netw.new <- subgraph.edges(graph=netw(object), 
                                              eids=E(netw(object))[preserved.edges],
                                              delete.vertices = FALSE)
              } else {
                netw.new<-object@netw
              }
              
              if(length(object@comm)!=0){
                comm.new <- object@comm
                comm.new$membership[!IDX] <- 0
                comm.new$modularity <- NA
              } else {
                comm.new <- object@comm
              }
              
              return(mgnet(abun=abun.new,
                           rela=rela.new,
                           norm=norm.new,
                           sample=object@sample,
                           taxa=object@taxa,
                           netw=netw.new,
                           comm=comm.new))
              # aggregate
              #-----------------------------------------#  
            } else if (trim=="aggregate") {
              
              # abun
              if(length(object@abun)!=0){
                abun.new<-object@abun[,IDX,drop=F]
                if(aggregate_to%in%colnames(abun.new)){
                  abun.new[,aggregate_to] <- abun.new[,aggregate_to] + rowSums(object@abun[,!IDX,drop=F])
                } else {
                  abun.new <- cbind(abun.new, aggregate_to=rowSums(object@abun[,!IDX,drop=F]))
                }
              } else {
                abun.new<-object@abun
              }
              # rela
              if(length(object@rela)!=0){
                rela.new<-object@rela[,IDX,drop=F]
                if(aggregate_to%in%colnames(rela.new)){
                  rela.new[,aggregate_to] <- rela.new[,aggregate_to] + rowSums(object@rela[,!IDX,drop=F])
                } else {
                  rela.new <- cbind(rela.new, aggregate_to=rowSums(object@rela[,!IDX,drop=F]))
                }
              } else {
                rela.new <- object@rela
              }
              # norm
              if(length(object@norm)!=0){
                norm.new<-object@norm[,IDX,drop=F]
                if(aggregate_to%in%colnames(norm.new)){
                  norm.new[,aggregate_to] <- norm.new[,aggregate_to] + rowSums(object@norm[,!IDX,drop=F])
                } else {
                  norm.new <- cbind(norm.new, aggregate_to=rowSums(object@norm[,!IDX,drop=F]))
                }
              } else {
                norm.new<-object@norm
              }
              # taxa
              if(length(object@taxa)!=0){
                taxa.new<-object@taxa[IDX,,drop=F]
                if(!(aggregate_to%in%rownames(taxa.new))){
                  taxa.new <- rbind(taxa.new,aggregate_to=rep(aggregate_to,ncol(object@taxa)))
                }
              } else {
                taxa.new<-object@taxa
              }
              # netw
              if(length(object@netw)!=0){
                netw.new<-igraph::subgraph(object@netw,IDX)
                if(!(aggregate_to%in%taxa_id(object))){
                  netw.new <- igraph::add_vertices(netw.new, nv=1, attr=list("name"=aggregate_to))
                }
              } else {
                netw.new<-object@netw
              }
              # comm
              if(length(object@comm)!=0){
                comm.new <- object@comm
                if(is.character(IDX)) IDX <- which(taxa_id(object)%in%IDX)
                comm.new$membership <- object@comm$membership[IDX]
                if(!(aggregate_to%in%taxa_id(object))){
                  comm.new$membership <- c(comm.new$membership,aggregate_to=0)
                }
                comm.new$vcount <- length(comm.new$membership)
                comm.new$modularity <- NA
              } else {
                comm.new <- object@comm
              }
              
              return(mgnet(abun=abun.new,
                           rela=rela.new,
                           norm=norm.new,
                           sample=object@sample,
                           taxa=taxa.new,
                           netw=netw.new,
                           comm=comm.new))
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