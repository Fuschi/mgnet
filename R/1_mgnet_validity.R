#------------------------------------------------------------------------------#
#' Internal: Assert Unique Row and Column Names
#'
#' Validates the presence and uniqueness of row and column names for a given object, 
#' either a matrix or a data.frame. This function is essential for maintaining data 
#' integrity within the `mgnet` class by ensuring that identifiers are both provided 
#' and uniquely identify each row and column. Error messages from the validation process 
#' are accumulated in an external vector, allowing for comprehensive feedback on data issues.
#'
#' @param obj The object to check, which must be either a matrix or a data.frame. 
#'        The function verifies that this object has non-null, unique row and column names.
#' @param errors A character vector accumulating error messages from validation checks.
#'
#' @return Returns the updated `errors` character vector. If no new validation issues are identified, 
#'         the vector is returned unchanged. If issues are found, corresponding error messages are 
#'         appended to the vector before it is returned, thereby accumulating messages across
#'         different validation checks.
#'
#' @keywords internal
.assertUniqueRowColNames <- function(obj, errors) {

  err <- character() 
  
  if (is.null(rownames(obj))) {
    err <- c(err, "row names are missing.")
  } else if (anyDuplicated(rownames(obj)) > 0) {
    err <- c(err, "row names must be unique.")
  }
  
  if (is.null(colnames(obj))) {
    err <- c(err, "column names are missing.")
  } else if (anyDuplicated(colnames(obj)) > 0) {
    err <- c(err, "column names must be unique.")
  }
  
  errors <- c(errors, err)
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Check for Reserved Keywords in Column Names
#'
#' This function checks if any of the specified reserved keywords are present as column names in a given object.
#' If such keywords are found, it appends an error message to the provided list of errors.
#' This check helps prevent conflicts with internal functionality that relies on these reserved keywords.
#'
#' @param obj The data object to check, expected to have column names (e.g., a data frame or a matrix with column names).
#' @param errors A character vector where error messages will be accumulated.
#'               New error messages are appended if reserved keywords are found as column names.
#'
#' @return A character vector of accumulated error messages, including any new errors found during this check.
#'
#' @keywords internal
.assertNoReservedKeywords <- function(obj, errors) {
  # Define the list of reserved keywords
  reservedKeywords <- c("taxa_id", "sample_sum")
  
  # Check if obj has column names to avoid unnecessary warnings
  if (!is.null(colnames(obj))) {
    # Iterate through each reserved keyword
    for (keyword in reservedKeywords) {
      if (keyword %in% colnames(obj)) {
        message <- switch(keyword,
                          "sample_id" = "column name 'sample_id' is a reserved keyword and cannot be used. This information could be encoded in the dimension names. See ?mgnet-class for details.",
                          "taxa_id" = "column name 'taxa_id' is a reserved keyword and cannot be used. This information could be encoded in the dimension names. See ?mgnet-class for details.",
        )
        errors <- c(errors, message)
      }
    }
  }
  
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Numeric Matrix
#'
#' Validates that the given object is a numeric matrix. This function is used 
#' internally to ensure data slots intended to be numeric matrices, such as abundance 
#' matrices, meet the required data type specifications.
#'
#' @param obj The object to check.
#' @param errors A character vector accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if `obj`
#'         is not a numeric matrix; otherwise, the original `errors` vector is returned.
#'
#' @keywords internal
.assertNumericMatrix <- function(obj, errors) {
  if (!is.matrix(obj) || !is.numeric(obj)) {
    errors <- c(errors, "must be a numeric matrix.")
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Data Frame
#'
#' Validates that the given object is a data frame. This function is essential for 
#' ensuring that slots expected to contain metadata or other tabular data are correctly 
#' formatted as data frames.
#'
#' @param obj The object to check.
#' @param errors A character vector for accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if `obj`
#'         is not a data frame; otherwise, the original `errors` vector is returned.
#'
#' @keywords internal
.assertDataFrame <- function(obj, errors) {
  if (!is.data.frame(obj)) {
    errors <- c(errors, "must be a data frame.")
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Character Matrix
#'
#' Validates that a given object is a character matrix. This function ensures that 
#' slots intended to hold character data, such as taxonomic classifications, adhere 
#' to the expected matrix format with character data types.
#'
#' @param obj The object to check.
#' @param errors A character vector for accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if `obj`
#'         is not a character matrix; otherwise, the original `errors` vector is returned.
#'
#' @keywords internal
.assertCharacterMatrix <- function(obj, errors) {
  if (!is.matrix(obj) || !all(apply(obj, c(1,2), is.character))) {
    errors <- c(errors, "must be a character matrix.")
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert All Elements Greater Than Zero
#'
#' Validates that all elements of a given numeric matrix or vector are greater than zero.
#' This function is crucial for slots where negative values are not permissible, ensuring
#' data integrity and consistency.
#'
#' @param obj The numeric matrix or vector to check.
#' @param errors A character vector for accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if any element
#'         in `obj` is not greater than zero; otherwise, the original `errors` vector is returned.
#'
#' @keywords internal
.assertAllPositive <- function(obj, errors) {
  if (any(obj < 0, na.rm = TRUE)) {
    errors <- c(errors, "all elements must be >= 0.")
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert igraph Network with Named Vertices
#'
#' Validates that the given object is an igraph network and checks that all vertices have names.
#' This function is crucial for ensuring network data integrity, requiring both the correct
#' structure and identifiable vertices for further analysis.
#'
#' @param obj The object to check, expected to be an igraph network.
#' @param errors A character vector for accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if `obj`
#'         is not an igraph network or if any vertices are unnamed; otherwise,
#'         the original `errors` vector is returned.
#'
#' @importFrom igraph V is_named
#' @keywords internal
.assertNamedIgraph <- function(obj, errors) {
  if (!inherits(obj, "igraph")) {
    errors <- c(errors, "must be an igraph network.")
  } else {
    if (any(is.na(V(obj)$name)) || !is_named(obj)) {
      errors <- c(errors, "all vertices in the igraph network must have names.")}
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert community Class Membership
#'
#' Validates that the given object is an instance of the `communities` class.
#' This function ensures that community detection results are stored in an appropriate
#' format, facilitating the analysis of microbial community structure.
#'
#' @param obj The object to check.
#' @param errors A character vector for accumulating error messages from validation checks.
#'
#' @return An updated character vector containing any new error messages if `obj`
#'         does not belong to the `communities` class; otherwise, the original
#'         `errors` vector is returned.
#'
#' @keywords internal
.assertCommunitiesClass <- function(obj, errors) {
  if (!inherits(obj, "communities")) {
    errors <- c(errors, "must belong to the communities class.")
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Column Names Between Two Slots
#'
#' Validates that two specified slots within an S4 object have the same number of columns
#' and matching column names. This function is crucial for ensuring consistency and
#' alignment between related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot to compare.
#' @param slotName2 The name of the second slot to compare.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If the number of columns or
#'         the column names between the two slots do not match, appropriate error messages
#'         are added to the vector. If there are no discrepancies, the original vector is
#'         returned unchanged.
#'
#' @keywords internal
.assertMatchingColumnNames <- function(obj, slotName1, slotName2, errors) {
  
  obj1 <- slot(obj, slotName1)
  obj2 <- slot(obj, slotName2)
  
  if ( ncol(obj1) != ncol(obj2) ){
    errors <- c(errors, sprintf("columns number of %s and %s are not equal.",
                                slotName1,slotName2))
  }else if ( any(colnames(obj1) != colnames(obj2)) ) {
    errors <- c(errors, sprintf("columns names of %s and %s do not match.",
                                slotName1, slotName2))
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Row Names Between Two Slots
#'
#' Validates that two specified slots within an S4 object have the same number of rows
#' and matching row names. This function is crucial for ensuring consistency and
#' alignment between related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot to compare.
#' @param slotName2 The name of the second slot to compare.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If the number of columns or
#'         the column names between the two slots do not match, appropriate error messages
#'         are added to the vector. If there are no discrepancies, the original vector is
#'         returned unchanged.
#'
#' @keywords internal
.assertMatchingRowNames <- function(obj, slotName1, slotName2, errors) {
  
  obj1 <- slot(obj, slotName1)
  obj2 <- slot(obj, slotName2)
  
  if ( nrow(obj1) != nrow(obj2) ){
    errors <- c(errors, sprintf("rows number of %s and %s are not equal.",
                                slotName1,slotName2))
  }else if ( any(rownames(obj1) != rownames(obj2)) ) {
    errors <- c(errors, sprintf("rows names of %s and %s do not match.",
                                slotName1, slotName2))
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Row and Column Names Between Two Slots
#'
#' Validates that two specified slots within an S4 object have the same number of 
#' rows and columns and matching row and column names. 
#' This function is crucial for ensuring consistency and
#' alignment between related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot to compare.
#' @param slotName2 The name of the second slot to compare.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If the number of columns or
#'         the column names between the two slots do not match, appropriate error messages
#'         are added to the vector. If there are no discrepancies, the original vector is
#'         returned unchanged.
#'
#' @keywords internal
.assertMatchingNames <- function(obj, slotName1, slotName2, errors) {
  
  obj1 <- slot(obj, slotName1)
  obj2 <- slot(obj, slotName2)
  
  if ( any(dim(obj1) != dim(obj2)) ){
    errors <- c(errors, sprintf("dimensions of %s and %s are not equal.",
                                slotName1,slotName2))
  } else if (length(errors)==0){
    if ( any(rownames(obj1) != rownames(obj2)) ) {
      errors <- c(errors, sprintf("rows names of %s and %s do not match.",
                                  slotName1, slotName2))}
    if ( any(colnames(obj1) != colnames(obj2)) ) {
      errors <- c(errors, sprintf("columns names of %s and %s do not match.",
                                  slotName1, slotName2))}
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Row Names and Column Names Between Two Slots
#' 
#' Validates that two specified slots within an S4 object have the first the same 
#' number of rows of the second columns and matching the respective names. 
#' This function is crucial for ensuring consistency and alignment between 
#' related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot to compare.
#' @param slotName2 The name of the second slot to compare.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If the number of columns or
#'         the column names between the two slots do not match, appropriate error messages
#'         are added to the vector. If there are no discrepancies, the original vector is
#'         returned unchanged.
#'
#' @keywords internal
.assertMatchingRowColsNames <- function(obj, slotName1, slotName2, errors) {
  
  obj1 <- slot(obj, slotName1)
  obj2 <- slot(obj, slotName2)
  
  if ( nrow(obj1) != ncol(obj2) ){
    errors <- c(errors, sprintf("rows number of %s and columns number %s are not equal.",
                                slotName1,slotName2))
  } else if ( any(rownames(obj1) != colnames(obj2)) ) {
    errors <- c(errors, sprintf("rows names of %s and columns names %s do not match.",
                                slotName1, slotName2))
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Col Names and Row Names Between Two Slots
#' 
#' Validates that two specified slots within an S4 object have the first the same 
#' number of columns of the second rows and matching the respective names. 
#' This function is crucial for ensuring consistency and alignment between 
#' related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot to compare.
#' @param slotName2 The name of the second slot to compare.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If the number of columns or
#'         the column names between the two slots do not match, appropriate error messages
#'         are added to the vector. If there are no discrepancies, the original vector is
#'         returned unchanged.
#'
#' @importFrom methods slot
#' @keywords internal
.assertMatchingColsRowsNames <- function(obj, slotName1, slotName2, errors) {
  
  obj1 <- slot(obj, slotName1)
  obj2 <- slot(obj, slotName2)
  
  if ( ncol(obj1) != nrow(obj2) ){
    errors <- c(errors, sprintf("columns number of %s and rows number %s are not equal.",
                                slotName1,slotName2))
  } else if ( any(colnames(obj1) != rownames(obj2)) ) {
    errors <- c(errors, sprintf("columns names of %s and rows names %s do not match.",
                                slotName1, slotName2))
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Dimension and Network Vertices Names
#' 
#' Validates that the specified slots within an S4 object match the names of the
#' network slot.
#' This function is crucial for ensuring consistency and alignment between 
#' related datasets within the `mgnet` class, such as ensuring that
#' matrices representing different types of data (e.g., raw and processed data) are
#' compatible in terms of their structure and labeling.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slot The name of the slot to compare.
#' @param which rows or columns of the slot to compare with the network.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. 
#'
#' @importFrom igraph V vcount
#' @importFrom methods slot
#' @keywords internal
.assertMatchingNamesVertices <- function(obj, slotName, which, errors) {
  
  netw <- slot(obj, "network")
  obj <- slot(obj, slotName)
  
  if(which=="rows"){
    slot_length <- nrow(obj)
    slot_names <- rownames(obj)
  } else if (which=="columns"){
    slot_length <- ncol(obj)
    slot_names <- colnames(obj)
  }
  
  netw_length <- vcount(netw)
  netw_names <- V(netw)$name
  
  if ( slot_length != netw_length ){
    errors <- c(errors, sprintf("%s number of %s and network vertices number are not equal.",
                                which,slotName))
  } else if ( any(slot_names != netw_names) ) {
    errors <- c(errors, sprintf("%s names of %s and network vertices names do not match.",
                                which, slotName))
  }
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Consistency Between Network and Community
#'
#' Validates the consistency between the network and its associated community within an S4 object. 
#' This function checks that a network is present and that the number of vertices in the network 
#' matches the number of community assignments. This is essential for ensuring that community 
#' detection results accurately reflect the structure of the network.
#'
#' @param obj An S4 object containing `network` and `community` slots.
#' @param errors A character vector that accumulates error messages from various validation checks. 
#'        New errors identified by this function are appended to this vector.
#'
#' @return Returns the updated `errors` character vector with any new error messages appended.
#'         If the network is missing or if the number of vertices in the network does not match 
#'         the number of community assignments, corresponding error messages are added.
#'
#' @details The function performs two main checks: first, it verifies that the `network` slot 
#'          is not empty, indicating that a network structure exists. Second, it compares the 
#'          number of vertices in the network (using `vcount`) with the length of the community 
#'          assignments in the `community$membership` vector. These validations ensure that 
#'          community data is aligned with the network's topology.
#'
#' @importFrom igraph vcount
#' @importFrom methods slot
#' @keywords internal
.assertMatchingCommunitiesNetwork <- function(obj, errors) {
  
  # Correcting 'object' to 'obj' based on the parameter name
  if(length(slot(obj, "network")) == 0){
    errors <- c(errors, "community cannot exist without the associated network.")
  } else if (length(slot(obj, "community")$membership) != vcount(slot(obj, "network"))){
    errors <- c(errors, "network and community slots must have the same number of vertices.")
  }
  
  return(errors)
}

#------------------------------------------------------------------------------#
# SET VALIDITY FUNCTION
#------------------------------------------------------------------------------#
# Set validity method for the mgnet class
# This method ensures that all instances of the mgnet class meet the predefined criteria
# for data integrity and consistency. Each slot of the mgnet object is checked for specific
# conditions, such as being non-empty, having the correct data type, and meeting domain-specific
# requirements (e.g., numeric matrices must have all elements >= 0).
setValidity("mgnet", function(object) {
  
  errors <- list()
  
  #CHECK ABUNDANCE
  #-------------------------------------#
  errors$abundance <- character()
  if( length(object@abundance)!=0 ){
    errors$abundance <- .assertNumericMatrix(object@abundance, errors$abundance)
    
    if ( length(errors$abundance)==0 ){
      errors$abundance <- .assertUniqueRowColNames(object@abundance, errors$abundance)
      errors$abundance <- .assertAllPositive(object@abundance, errors$abundance)
      errors$abundance <- .assertNoReservedKeywords(object@abundance, errors$abundance)
    }
  }
  
  #CHECK LOG-ABUNDANCE
  #-------------------------------------#
  errors$log_abundance <- character()
  if( length(object@log_abundance)!=0 ){
    errors$log_abundance <- .assertNumericMatrix(object@log_abundance, errors$log_abundance)
    
    if ( length(errors$log_abundance)==0 ){
      errors$log_abundance <- .assertUniqueRowColNames(object@log_abundance, errors$log_abundance)
      errors$log_abundance <- .assertNoReservedKeywords(object@log_abundance, errors$log_abundance)
      
    }
  }
      
  #CHECK INFO_SAMPLE
  #-------------------------------------#
  errors$info_sample <- character()
  if( length(object@info_sample)!=0 ){
    errors$info_sample <- .assertDataFrame(object@info_sample, errors$info_sample)
    
    if ( length(errors$info_sample)==0 ){
      errors$info_sample <- .assertUniqueRowColNames(object@info_sample, errors$info_sample)
      errors$info_sample <- .assertNoReservedKeywords(object@info_sample, errors$info_sample)
    }
  } 
  
  #CHECK LINEAGE
  #-------------------------------------#
  errors$lineage <- character()
  if( length(object@lineage)!=0 ){
    errors$lineage <- .assertCharacterMatrix(object@lineage, errors$lineage)
    
    if ( length(errors$lineage)==0 ){
      errors$lineage <- .assertUniqueRowColNames(object@lineage, errors$lineage)
      errors$lineage <- .assertNoReservedKeywords(object@lineage, errors$lineage)
    }
  }
  
  #CHECK INFO_TAXA
  #-------------------------------------#
  errors$info_taxa <- character()
  if( length(object@info_taxa)!=0 ){
    errors$info_taxa <- .assertDataFrame(object@info_taxa, errors$info_taxa)
    
    if ( length(errors$info_taxa)==0 ){
      errors$info_taxa <- .assertUniqueRowColNames(object@info_taxa, errors$info_taxa)
      errors$info_taxa <- .assertNoReservedKeywords(object@info_taxa, errors$info_taxa)
    }
  }
  
  # CHECK NETWORK
  #-------------------------------------#
  errors$network <- character()
  if(length(object@network)!=0){
    errors$network <- .assertNamedIgraph(object@network, errors$network)
  }
  
  # CHECK community
  #-------------------------------------#
  errors$community <- character()
  if(length(object@community)!=0){
    errors$community <- .assertNamedIgraph(object@community, errors$community)
  }
  
  # CHECK RECIPROCAL PROPERTIES
  #-------------------------------------#
  errors$reciprocal <- character()
  
  if(length(object@abundance)!=0 && length(errors$abundance)==0){
    if(length(object@log_abundance)!=0 && length(errors$log_abundance)==0){
      errors$reciprocal <- .assertMatchingNames(object, "abundance", "log_abundance", errors$reciprocal)}
    if(length(object@info_sample)!=0 && length(errors$info_sample)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "abundance", "info_sample", errors$reciprocal)}
    if(length(object@lineage)!=0 && length(errors$lineage)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "abundance", "lineage", errors$reciprocal)}
    if(length(object@info_taxa)!=0 && length(errors$info_taxa)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "abundance", "info_taxa", errors$reciprocal)}
    if(length(object@network)!=0 && length(errors$network)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "abundance", "rows", errors$reciprocal)}
  }
  
  if(length(object@log_abundance)!=0 && length(errors$log_abundance)==0){
    if(length(object@info_sample)!=0 && length(errors$info_sample)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "log_abundance", "info_sample", errors$reciprocal)}
    if(length(object@lineage)!=0 && length(errors$lineage)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "log_abundance", "lineage", errors$reciprocal)}
    if(length(object@info_taxa)!=0 && length(errors$info_taxa)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "log_abundance", "info_taxa", errors$reciprocal)}
    if(length(object@network)!=0 && length(errors$network)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "log_abundance", "rows", errors$reciprocal)}
  }
  
  if(length(object@lineage)!=0 && length(errors$lineage)==0){
    if(length(object@info_taxa)!=0 && length(errors$info_taxa)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "lineage", "info_taxa", errors$reciprocal)}
    if(length(object@network)!=0 && length(errors$network)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "lineage", "rows", errors$reciprocal)}
  }
  
  if(length(object@info_taxa)!=0 && length(errors$info_taxa)==0){
    if(length(object@network)!=0 && length(errors$network)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "info_taxa", "rows", errors$reciprocal)}
  }
  
  if(length(object@community)!=0 && length(errors$community)==0){
    errors$reciprocal <- .assertMatchingCommunitiesNetwork(object, errors$reciprocal)
  }

  # FORMAT ERROR OUTPUT
  errors <- Filter(function(x) length(x) > 0, errors)
  if(length(errors)==0){
    return(TRUE)
  } else {
    errors <- lapply(errors, function(x) paste("- ", x, "\n", sep="") )
    errors <- lapply(errors, function(x) paste(x, collapse=""))
    errors <- mapply(function(slot,msg){
      paste("\nDEBUGGER of ",slot,":\n", msg, sep="")
    }, slot=names(errors), msg=errors)
    errors <- paste(c(errors), sep="\n")
    errors <- gsub("reciprocal", "reciprocal properties", errors)
    return(errors)
  }
})
#------------------------------------------------------------------------------#
# END SET VALIDITY FUNCTION
#------------------------------------------------------------------------------#