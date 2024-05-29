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
  
  if(length(object@info_sample)==0){
    return(object[numeric(0),])
  }
  
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Filter the info_sample data frame based on the conditions
  filtered_info_sample <- dplyr::filter(info_sample(object, .fmt = "tbl"), !!!conditions)
  
  # Check if the filtered result is empty
  if(nrow(filtered_info_sample) == 0) {
    return(object[integer(0),])
  }
  
  # Extract the sample_id of the filtered samples
  filtered_sample_ids <- filtered_info_sample$sample_id
  
  # Subset the object using the mgnet extractor method based on filtered sample IDs
  subsetted_object <- object[filtered_sample_ids, ]
  
  return(subsetted_object)
})


setMethod("filter_info_sample", "mgnetList", function(object, ...) {
  # Capture the filtering conditions as quosures
  conditions <- rlang::enquos(...)
  
  # Apply the filter_samples method to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    
    if(length(mgnet_obj@info_sample)==0){
      return(mgnet_obj[numeric(0),])
    }
    
    # Filter the info_sample data frame based on the conditions
    filtered_info_sample <- dplyr::filter(info_sample(mgnet_obj, .fmt = "tbl"), !!!conditions)
    
    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_info_sample) == 0) {
      return(object[integer(0),])
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
#' merged dataset from `info_taxa`, `lineage` and the communities membership from `community` slots.
#' @param trim How to handle taxa not meeting criteria: "yes" to remove them,
#'        "no" to set their abundance to zero (min sample value for norm_abundance), 
#'        or "aggregate" to combine them into a single category.
#' @param aggregate_to Name for the aggregated taxa category when using `trim = "aggregate"`.
#'        Defaults to "aggregate". 
#'
#' @return An `mgnet` or `mgnetList` object with the taxa subset according to the
#' specified conditions. The structure of the object remains intact, but taxa that
#' do not meet the conditions are removed from the abundance data, lineage information,
#' and potentially network and community analysis results.
#'
#' @details The function first merges `info_taxa` data with `lineage` data and 
#' communities membership from `community`to create a comprehensive taxonomic dataset. 
#' It then filters this dataset based on the provided conditions, identifying which taxa to retain. 
#' The resulting filtered taxa list is used to subset the original `mgnet` or `mgnetList` object, 
#' thus affecting related abundance and network data but not altering the sample metadata.
#'
#' @export
#' @importFrom dplyr filter inner_join
#' @importFrom methods validObject
#' @importFrom rlang enquos
#' @name filter_info_taxa
#' @aliases filter_info_taxa,mgnet-method filter_info_taxa,mgnetList-method
#' @seealso \code{\link[dplyr]{filter}} for details on filter expressions.
setGeneric("filter_info_taxa", function(object, ...,
                                        trim = "yes", aggregate_to = NULL) {standardGeneric("filter_info_taxa")})

setMethod("filter_info_taxa", "mgnet", function(object, ..., trim = "yes", aggregate_to = NULL) {
  
  if(length(object@info_taxa)==0 & length(object@lineage)==0 && length(object@community)==0){
    return(object[,numeric(0)])
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
  
  if(length(object@info_taxa)==0 & length(object@lineage)==0 & length(community(object)) !=0 ) {
    
    merged_taxa <- tibble("taxa_id" = taxa_id(object), "community" = community_members(object))
    
  } else if(length(community(object)) !=0 ){
    
    merged_taxa <- merged_taxa %>%
      left_join(tibble("taxa_id" = taxa_id(object), "community" = community_members(object)),
                by = join_by(taxa_id))
  }
  
  
  # Apply conditions to filter taxa
  filtered_taxa <- dplyr::filter(merged_taxa, !!!conditions)
  
  # Check if the filtered result is empty
  if(nrow(filtered_taxa) == 0) {
    return(object[,integer(0)])
  }
  
  # Extract taxa IDs that meet the filtering conditions
  final_pass <- filtered_taxa$taxa_id
  
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

setMethod("filter_info_taxa", "mgnetList", function(object, ..., trim = "yes", aggregate_to = NULL) {
  
  # Apply filtering to each mgnet object within the mgnetList
  object@mgnets <- sapply(object@mgnets, function(mgnet_obj) {
    filtered_taxa_obj <- filter_info_taxa(object = mgnet_obj, ... = ..., 
                                          trim = trim, aggregate_to = aggregate_to)
    
    # Capture filtering conditions as quosures
    conditions <- rlang::enquos(...)
    
    # Check if the filtered result is empty and throw an error if so
    if(length(filtered_taxa_obj) == 0) { # Adjust this condition as needed based on actual implementation
      return(object[,integer(0)])
    }
    
    return(filtered_taxa_obj)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  validObject(object)
  return(object)
})