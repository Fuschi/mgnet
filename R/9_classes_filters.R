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
#' @param trim A parameter to handle taxa not meeting criteria: "yes" to remove them,
#'        "no" to set their abundance to zero, or "aggregate" to combine them into a 'merged' category.
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
setGeneric("refine_taxa", function(object, ..., condition = "AND", trim = "yes") standardGeneric("refine_taxa"))

setMethod("refine_taxa", "mgnet",
          function(object,...,condition="AND",trim="yes"){
            
            condition <- match.arg(condition, c("AND","OR"))
            trim <- match.arg(trim, c("yes","no","aggregate"))
            
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
              
              ifelse(length(object@abundance)!=0, abundance.new<-object@abundance[,IDX,drop=F], abundance.new<-object@abundance)
              ifelse(length(object@lineage)!=0, lineage.new<-object@lineage[IDX,,drop=F], lineage.new<-object@lineage)
              ifelse(length(object@info_taxa)!=0, info_lineage.new<-object@info_taxa[IDX,,drop=F], info_lineage.new<-object@info_taxa)
              ifelse(length(object@log_abundance)!=0, log_abundance.new<-object@log_abundance[,IDX,drop=F], log_abundance.new<-object@log_abundance)
              
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,IDX)
              } else {
                network.new<-object@network
              }
              
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(IDX)) IDX <- which(taxa_id(object)%in%IDX)
                community.new$membership <- object@community$membership[IDX]
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
              
              # no
              #-----------------------------------------#
            } else if (trim=="no") {
              
              if(length(object@abundance)!=0){
                abundance.new <- object@abundance
                abundance.new[,!IDX] <- 0
              } else {
                abundance.new <- object@abundance
              }
              
              if(length(object@log_abundance)!=0){
                log_abundance.new <- object@log_abundance
                sample_min <- apply(object@log_abundance, 1, min)
                log_abundance.new[,!IDX] <- sample_min # This could be an error (see later)
              } else {
                log_abundance.new <- object@log_abundance
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
                abundance.new<-object@abundance[,IDX,drop=F]
                if("merged"%in%colnames(abundance.new)){
                  abundance.new[,"merged"] <- abundance.new[,"merged"] + rowSums(object@abundance[,!IDX,drop=F])
                } else {
                  abundance.new <- cbind(abundance.new, "merged"=rowSums(object@abundance[,!IDX,drop=F]))
                }
              } else {
                abundance.new<-object@abundance
              }
              # taxa
              if(length(object@lineage)!=0){
                lineage.new<-object@lineage[IDX,,drop=F]
                if(!("merged"%in%rownames(lineage.new))){
                  lineage.new <- rbind(lineage.new,"merged"=rep("merged",length(ranks(object))))
                }
              } else {
                lineage.new<-object@lineage
              }
              # info_taxa
              if(length(object@info_taxa)!=0){
                info_lineage.new<-object@info_taxa[IDX,,drop=F]
                if(!("merged"%in%rownames(info_lineage.new))){
                  info_lineage.new <- rbind(info_lineage.new,"merged"=rep("merged",ncol(object@info_taxa)))
                }
              } else {
                info_lineage.new<-object@info_taxa
              }
              # log_data
              if(length(object@log_abundance)!=0){
                log_abundance.new<-object@log_abundance[,IDX,drop=F]
                if("merged"%in%colnames(log_abundance.new)){
                  log_abundance.new[,"merged"] <- log_abundance.new[,"merged"] + rowSums(object@log_abundance[,!IDX,drop=F])
                } else {
                  log_abundance.new <- cbind(log_abundance.new, "merged"=rowSums(object@log_abundance[,!IDX,drop=F]))
                }
              } else {
                log_abundance.new<-object@log_abundance
              }
              # netw
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,IDX)
                if(!("merged"%in%taxa_id(object))){
                  network.new <- igraph::add_vertices(network.new, nv=1, attr=list("name"="merged"))
                }
              } else {
                network.new<-object@network
              }
              # comm
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(IDX)) IDX <- which(taxa_id(object)%in%IDX)
                community.new$membership <- object@community$membership[IDX]
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
            
          })


setMethod("refine_taxa", "mgnetList",
          function(object, ..., condition="AND", trim="yes") {
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


# FILTER INFO SAMPLE
#------------------------------------------------------------------------------#
#' Filter Samples in mgnet Objects Based on info_sample Conditions
#'
#' @description
#' This function filters samples in an `mgnet` object or each `mgnet` object within
#' an `mgnetList` based on specified conditions applied to the `info_sample` data frame.
#' The conditions for filtering are based on expressions applied to columns within
#' the `info_sample` data frame.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Conditions to filter samples, which are passed unquoted and can
#' include dplyr-style filtering expressions. These conditions are applied to the
#' columns of the `info_sample` data frame, which is accessed as a tibble with
#' `sample_id` serving as the key column.
#'
#' @return An `mgnet` or `mgnetList` object with samples filtered according to the
#' specified conditions. The filtered object(s) retain only those samples that meet
#' the filtering criteria specified by the `...` arguments.
#'
#' @details The `filter_info_sample` function leverages the `info_sample` getter function
#' to access the metadata of samples as a tibble, where the row names are transformed
#' into a `sample_id` column. This enables seamless integration with dplyr's `filter`
#' function for specifying complex filtering conditions. The function then uses the
#' `sample_id` column to identify which samples should be retained in the filtered `mgnet`
#' or `mgnetList` object.
#'
#' @export
#' @importFrom dplyr filter
#' @importFrom rlang enquos
#' @name filter_info_sample
#' @aliases filter_info_sample,mgnet-method filter_info_sample,mgnetList-method
#' @seealso \code{\link[dplyr]{filter}} for details on filter conditions.
setGeneric("filter_info_sample", function(object, ...) {
  standardGeneric("filter_info_sample")
})

setMethod("filter_info_sample", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)

  # Filter the info_sample data frame based on the conditions
  filtered_info_sample <- dplyr::filter(info_sample(object, .fmt = "tbl"), !!!conditions)

  # Check if the filtered result is empty
  if(nrow(filtered_info_sample) == 0) {
    stop("No samples meet the specified filtering conditions. Please revise your criteria.")
  }

  # Extract the sample_id of the filtered samples
  filtered_sample_ids <- filtered_info_sample$sample_id

  # Subset the object using the mgnet extractor method based on filtered sample IDs
  subsetted_object <- object[filtered_sample_ids,]

  return(subsetted_object)
})


setMethod("filter_info_sample", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)

  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {

    # Filter the info_sample data frame based on the conditions
    filtered_info_sample <- dplyr::filter(info_sample(mgnet_obj, .fmt = "tbl"), !!!conditions)

    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_info_sample) == 0) {
      stop("No samples meet the specified filtering conditions in one or more mgnet objects. Please revise your criteria.")
    }

    # Extract the sample_id of the filtered samples
    filtered_sample_ids <- filtered_info_sample$sample_id

    # Subset the object using the mgnet extractor method based on filtered sample IDs
    subsetted_object <- mgnet_obj[filtered_sample_ids,]

    return(subsetted_object)

  }, simplify = FALSE, USE.NAMES = TRUE)

  return(object)
})


# FILTER INFO TAXA
#------------------------------------------------------------------------------#
#' Filter Taxa in mgnet Objects Based on Taxonomic and Metadata Conditions
#'
#' @description
#' This function filters taxa in an `mgnet` object or each `mgnet` object within
#' an `mgnetList` based on specified conditions applied to the merged data from
#' the `info_taxa` and `lineage` slots. This enables filtering taxa using both
#' taxonomic classification and additional taxa metadata.
#'
#' @param object An `mgnet` or `mgnetList` object.
#' @param ... Conditions for filtering taxa, passed unquoted and allowing for
#' dplyr-style filtering expressions. These conditions apply to columns in the
#' merged dataset from `info_taxa` and `lineage`.
#'
#' @return An `mgnet` or `mgnetList` object with the taxa subset according to the
#' specified conditions. The structure of the object remains intact, but taxa that
#' do not meet the conditions are removed from the abundance data, lineage information,
#' and potentially network and community analysis results.
#'
#' @details The function first merges `info_taxa` data with `lineage` data to create
#' a comprehensive taxonomic dataset. It then filters this dataset based on the
#' provided conditions, identifying which taxa to retain. The resulting filtered
#' taxa list is used to subset the original `mgnet` or `mgnetList` object, thus
#' affecting related abundance and network data but not altering the sample metadata.
#'
#' @export
#' @importFrom dplyr filter inner_join
#' @importFrom methods validObject
#' @importFrom rlang enquos
#' @name filter_info_taxa
#' @aliases filter_info_taxa,mgnet-method filter_info_taxa,mgnetList-method
#' @seealso \code{\link[dplyr]{filter}} for details on filter expressions.
setGeneric("filter_info_taxa", function(object, ...) {standardGeneric("filter_info_taxa")})

setMethod("filter_info_taxa", "mgnet", function(object, ...) {

  if(length(object@info_taxa)==0 & length(object@lineage)==0){
    stop("slots info_taxa and lineage missing")
  }

  # Capture filtering conditions as quosures
  conditions <- rlang::enquos(...)

  # Merge 'info_taxa' and 'lineage' into a comprehensive dataset
  if(length(object@info_taxa)!=0 & length(object@lineage)!=0){

    merged_taxa <- inner_join(info_taxa(object, .fmt = "tbl"),
                             lineage(object, .fmt = "tbl"),
                             by = "taxa_id")

  } else if (length(object@info_taxa)==0 & length(object@lineage)!=0) {

    merged_taxa <- lineage(object, .fmt = "tbl")

  } else if (length(object@info_taxa)!=0 & length(object@lineage)==0) {

    merged_taxa <- info_taxa(object, .fmt = "tbl")

  }


  # Apply conditions to filter taxa
  filtered_taxa <- dplyr::filter(merged_taxa, !!!conditions)

  # Check if the filtered result is empty
  if(nrow(filtered_taxa) == 0) {
    stop("No taxa meet the specified filtering conditions. Please revise your criteria.")
  }

  # Extract taxa IDs that meet the filtering conditions
  filtered_taxa_ids <- filtered_taxa$taxa_id

  # Subset the mgnet object to retain only filtered taxa
  subsetted_object <- object[, filtered_taxa_ids]

  return(subsetted_object)
})

setMethod("filter_info_taxa", "mgnetList", function(object, ...) {
  # Capture filtering conditions as quosures for list application
  conditions <- rlang::enquos(...)

  # Apply filtering to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    filtered_taxa_obj <- filter_info_taxa(mgnet_obj, !!!conditions)

    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_taxa_obj) == 0) { # Adjust this condition as needed based on actual implementation
      stop("No taxa meet the specified filtering conditions in one or more mgnet objects. Please revise your criteria.")
    }

    return(filtered_taxa_obj)
  }, simplify = FALSE, USE.NAMES = TRUE)

  validObject(object)
  return(object)
})

 
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
            if(!is.null(abundance_criteria) && length(object@abundance)==0) stop("cannot set abundace_criteria without the abundance matrix")
            if(!is.null(relative_criteria) && !("sample_sum"%in%info_sample_vars(object))) {
              stop("The 'sample_sum' column is missing in 'info_sample'. ",
                   "Please ensure 'sample_sum' is calculated and included. ",
                   "This can be done automatically by creating or updating the mgnet object ",
                   "with valid abundance data, or through the 'update_sample_sum()' function.",
                   "\nSee '?update_sample_sum' for more information.")
            }
            if(!is.null(log_abundance_criteria) && length(object@log_abundance)==0) stop("cannot set log_abundace_criteria without the log_abundance matrix")
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(log_abundance_criteria)) stop("cannot missing all criteria")

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
            relative_pass <- if(length(relative(object))!=0){
              apply_criteria(relative(object), relative_criteria, condition)
            } else {
              rep(TRUE, nsample(object))
            }
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

            # Filter samples based on the combined criteria
            if (!any(final_pass)) {
              stop("No samples meet the specified criteria.")
            }

            subsetted_object <- object[final_pass, , drop = FALSE]
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
            if(!is.null(abundance_criteria) && length(object@abundance)==0) stop("cannot set abundace_criteria without the abundance matrix")
            if(!is.null(relative_criteria) && !("sample_sum"%in%info_sample_vars(object))) {
              stop("The 'sample_sum' column is missing in 'info_sample'. ",
                   "Please ensure 'sample_sum' is calculated and included. ",
                   "This can be done automatically by creating or updating the mgnet object ",
                   "with valid abundance data, or through the 'update_sample_sum()' function.",
                   "\nSee '?update_sample_sum' for more information.")
            }
            if(!is.null(log_abundance_criteria) && length(object@log_abundance)==0) stop("cannot set log_abundace_criteria without the log_abundance matrix")
            if(!is.null(abundance_criteria) & !is.null(relative_criteria) & !is.null(log_abundance_criteria)) stop("cannot missing all criteria")

            # Automatically remove taxa that are always equal to zero across all samples
            if(exclude_absent_taxa && length(abundance) != 0) {
              zero_taxa <- colSums(abundance(object)) != 0
              if(all(zero_taxa==0)) stop("all elements in abundance cannot be equal to zero")
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

            # Filter samples based on the combined criteria
            if (!any(final_pass)) {
              stop("No samples meet the specified criteria.")
            }
            
            
            # yes
            #-----------------------------------------#
            if(trim=="yes"){
              
              ifelse(length(object@abundance)!=0, abundance.new<-object@abundance[,final_pass,drop=F], abundance.new<-object@abundance)
              ifelse(length(object@lineage)!=0, lineage.new<-object@lineage[final_pass,,drop=F], lineage.new<-object@lineage)
              ifelse(length(object@info_taxa)!=0, info_lineage.new<-object@info_taxa[final_pass,,drop=F], info_lineage.new<-object@info_taxa)
              ifelse(length(object@log_abundance)!=0, log_abundance.new<-object@log_abundance[,final_pass,drop=F], log_abundance.new<-object@log_abundance)
              
              if(length(object@network)!=0){
                network.new<-igraph::subgraph(object@network,final_pass)
              } else {
                network.new<-object@network
              }
              
              if(length(object@community)!=0){
                community.new <- object@community
                if(is.character(final_pass)) final_pass <- which(taxa_id(object)%in%final_pass)
                community.new$membership <- object@community$membership[final_pass]
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


# SELECT INFO SAMPLE
#------------------------------------------------------------------------------#
#' Select Columns from info_sample in mgnet Objects
#'
#' This function allows users to select specific columns from the `info_sample` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's select semantics to provide a flexible interface for column selection based on column names or conditions.
#'
#' @description
#' `select_info_sample` uses dplyr's select functionality to enable precise selection of columns from the `info_sample`
#' data frame in `mgnet` objects. This can be particularly useful for simplifying metadata before further analysis
#' or for focusing on a subset of metadata attributes.
#'
#' @param object An `mgnet` or `mgnetList` object from which columns will be selected.
#' @param ... Conditions specifying which columns to select, passed to dplyr::select().
#'        This can include a variety of selectors like column names, dplyr helper functions (starts_with, ends_with, contains, etc.), or indices.
#'        For detailed usage, refer to \code{\link[dplyr]{select}}.
#'
#' @return An `mgnet` object with the `info_sample` slot updated to include only the selected columns, or an `mgnetList` object
#'         where each contained `mgnet` object has its `info_sample` slot similarly updated.
#'
#' @seealso
#' \code{\link[dplyr]{select}} for details on the selection conditions.
#' \code{\link{mgnet}} and \code{\link{mgnetList}} for details on the object structures.
#'
#' @export
#' @name select_info_sample
#' @aliases select_info_sample,mgnet-method select_info_sample,mgnetList-method
#' @importFrom dplyr select
#' @importFrom rlang enquos
setGeneric("select_info_sample", function(object, ...) {standardGeneric("select_info_sample")})

setMethod("select_info_sample", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  # Apply dplyr::select to info_sample
  selected_info_sample <- dplyr::select(info_sample(object, .fmt = "df"), !!!conditions)

  if(ncol(selected_info_sample) == 0) {
    stop("No columns meet the specified selecting conditions. Please revise your criteria.")
  }

  # Update info_sample directly
  object@info_sample <- selected_info_sample 
  validObject(object)
  
  return(object)
})


setMethod("select_info_sample", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)

  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {

    # Apply dplyr::select to info_sample
    selected_info_sample <- dplyr::select(info_sample(mgnet_obj, .fmt = "df"), !!!conditions)

    if(ncol(selected_info_sample) == 0) {
      stop("No columns meet the specified selecting conditions. Please revise your criteria.")
    }

    # Update info_sample directly
    mgnet_obj@info_sample <- selected_info_sample
    validObject(mgnet_obj)

    return(mgnet_obj)

  }, simplify = FALSE, USE.NAMES = TRUE)

  validObject(object)
  return(object)
})

# SELECT INFO TAXA
#------------------------------------------------------------------------------#
#' Select Columns from info_taxa in mgnet Objects
#'
#' This function allows users to select specific columns from the `info_taxa` slot of `mgnet` or `mgnetList` objects.
#' It utilizes dplyr's select semantics to provide a flexible interface for column selection based on column names or conditions.
#'
#' @description
#' `select_info_taxa` uses dplyr's select functionality to enable precise selection of columns from the `info_taxa`
#' data frame in `mgnet` objects. This can be particularly useful for simplifying metadata before further analysis
#' or for focusing on a subset of metadata attributes.
#'
#' @param object An `mgnet` or `mgnetList` object from which columns will be selected.
#' @param ... Conditions specifying which columns to select, passed to dplyr::select().
#'        This can include a variety of selectors like column names, dplyr helper functions (starts_with, ends_with, contains, etc.), or indices.
#'        For detailed usage, refer to \code{\link[dplyr]{select}}.
#'
#' @return An `mgnet` object with the `info_taxa` slot updated to include only the selected columns, or an `mgnetList` object
#'         where each contained `mgnet` object has its `info_taxa` slot similarly updated.
#'
#' @seealso
#' \code{\link[dplyr]{select}} for details on the selection conditions.
#' \code{\link{mgnet}} and \code{\link{mgnetList}} for details on the object structures.
#'
#' @export
#' @name select_info_taxa
#' @aliases select_info_taxa,mgnet-method select_info_taxa,mgnetList-method
#' @importFrom dplyr select
#' @importFrom rlang enquos
setGeneric("select_info_taxa", function(object, ...) {standardGeneric("select_info_taxa")})

setMethod("select_info_taxa", "mgnet", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  # Apply dplyr::select to info_taxa
  selected_info_taxa <- dplyr::select(info_taxa(object, .fmt = "df"), !!!conditions)
  
  if(ncol(selected_info_taxa) == 0) {
    stop("No columns meet the specified selecting conditions. Please revise your criteria.")
  }
  
  # Update info_sample directly
  object@info_taxa <- selected_info_taxa
  validObject(object)
  
  return(object)
})


setMethod("select_info_taxa", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    
    # Apply dplyr::select to info_sample
    selected_info_taxa <- dplyr::select(info_taxa(mgnet_obj, .fmt = "df"), !!!conditions)
    
    if(ncol(selected_info_taxa) == 0) {
      stop("No columns meet the specified selecting conditions. Please revise your criteria.")
    }
    
    # Update info_sample directly
    mgnet_obj@info_taxa <- selected_info_taxa
    validObject(mgnet_obj)
    
    return(mgnet_obj)
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  validObject(object)
  return(object)
})


# SPLIT MGNET
#------------------------------------------------------------------------------#
#' Split an mgnet Object into an mgnetList Based on Specified Columns
#'
#' @description
#' This method splits an `mgnet` object into subsets based on unique combinations
#' of values in specified columns from the `info_sample` data frame. Each subset
#' corresponds to an `mgnet` object in the resulting `mgnetList`, facilitating
#' analyses specific to each unique combination of metadata values.
#'
#' @param object An `mgnet` object.
#' @param ... Columns in the `info_sample` data frame used to define subsets.
#'        You can specify columns by name, using any of the selection methods supported by `dplyr::select()`.
#' @return An `mgnetList` object containing `mgnet` objects for each unique combination
#'         of specified column values.
#'
#' @details
#' The `split_mgnet` function leverages `dplyr` functionality to identify unique combinations
#' of specified column values in the `info_sample` data frame. It then creates a new `mgnet`
#' object for each combination, allowing for tailored analysis per group.
#'
#' This function is particularly useful in exploratory data analysis and preprocessing stages
#' where data need to be examined or analyzed based on specific grouping variables.
#'
#' @note
#' This function depends on `dplyr` for data manipulation. Ensure that `dplyr` is installed
#' and loaded in your R session.
#'
#' @export
#' @aliases split_mgnet,mgnet-method
#' @importFrom dplyr select distinct group_by group_split semi_join
#' @importFrom rlang enquos
setGeneric("split_mgnet", function(object, ...) {
  standardGeneric("split_mgnet")
})

setMethod("split_mgnet", "mgnet", function(object, ...) {
  # Capture the column names as quosures
  split_cols <- enquos(...)

  # Get the unique combinations from info_sample for the specified columns
  info_sample <- info_sample(object, "tbl")
  groups <- info_sample %>%
    dplyr::select(!!!split_cols) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!!split_cols) %>%
    dplyr::group_split()

  # Initialize an empty list to store the subsetted mgnet objects
  mgnet_list <- mgnetList()

  # For each group, filter the info_sample based on the unique combination and create a new mgnet object
  for (i in seq_along(groups)) {
    # Filter info_sample based on the group
    group_df <- groups[[i]]
    filtered_info_sample <- dplyr::semi_join(info_sample, group_df, by = names(group_df))

    # Filter the original mgnet object to create a subset based on the filtered info_sample
    subsetted_mgnet <- object[filtered_info_sample$sample_id, , drop = FALSE]

    # Add the subsetted mgnet object to the list with a meaningful name
    group_name <- paste(group_df %>% dplyr::pull(names(group_df)) %>% unlist(), collapse = "-")
    mgnet_list[[group_name]] <- subsetted_mgnet
  }

  return(mgnet_list)
})


#' # MGNET_LONGER
#' #------------------------------------------------------------------------------#
#' #' Convert mgnet Object to Longer Format
#' #'
#' #' @description 
#' #' Converts data from an `mgnet` object into a longer format suitable for easy analysis and casting,
#' #' combining abundance, sample metadata, taxonomic information, and other relevant data.
#' #'
#' #' @usage mgnet_longer(object)
#' #'
#' #' @param object An `mgnet` object.
#' #'
#' #' @return A `tibble` in a long format containing the combined data from the `mgnet` object.
#' #'
#' #' @export
#' #' @examples
#' #' # Assuming `mgnet_obj` is an existing `mgnet` object
#' #' longer_df <- mgnet_longer(mgnet_obj)
#' #'
#' #' @importFrom tidyr pivot_longer
#' #' @importFrom dplyr left_join
#' setGeneric("mgnet_longer", function(object, 
#'                                     abundance = FALSE, relative = FALSE, log_abundance = FALSE,
#'                                     info_sample = NULL, info_taxa = NULL, lineage = NULL) standardGeneric("mgnet_longer"))
#' 
#' setMethod("mgnet_longer", "mgnet", function(object){
#'   if(length(object@abundance) == 0) {
#'     stop("The 'abundance' slot must be present in the 'mgnet' object.")
#'   }
#'   
#'   # Convert abundance matrix to long format
#'   abundance_long <-
#'   
#'   # Join with sample metadata if available
#'   if(length(object@info_sample) != 0) {
#'     sample_info_long <- pivot_longer(as.data.frame(t(object@info_sample)), 
#'                                      cols = everything(), names_to = "sampleID", 
#'                                      values_to = "sample_info")
#'     abundance_long <- left_join(abundance_long, sample_info_long, by = "sampleID")
#'   }
#'   
#'   # Join with taxa metadata if available
#'   if(length(object@info_taxa) != 0) {
#'     taxa_info_long <- pivot_longer(as.data.frame(t(object@info_taxa)), 
#'                                    cols = everything(), names_to = "taxaID", 
#'                                    values_to = "taxa_info")
#'     abundance_long <- left_join(abundance_long, taxa_info_long, by = "taxaID")
#'   }
#'   
#'   # Add relative abundance if available
#'   if("sample_sum" %in% names(object@info_sample)) {
#'     abundance_long <- abundance_long %>%
#'       mutate(relative_abundance = abundance / object@info_sample$sample_sum[sampleID])
#'   }
#'   
#'   # Add log abundance if available
#'   if(length(object@log_abundance) != 0) {
#'     log_abundance_long <- as.data.frame(as.table(object@log_abundance)) %>%
#'       rename(sampleID = Var1, taxaID = Var2, log_abundance = Freq)
#'     abundance_long <- left_join(abundance_long, log_abundance_long, by = c("sampleID", "taxaID"))
#'   }
#'   
#'   # Return as a tibble for consistency with tidyverse output
#'   return(as_tibble(abundance_long))
#' })
