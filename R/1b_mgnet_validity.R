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
#' Internal: Check for Reserved Keywords in Column Names or Vertex Names of a Specified Slot
#'
#' This function checks if any of the specified reserved keywords are present as column names in a matrix or data frame,
#' or as vertex names in an igraph object within a specified slot of an S4 object.
#' If such keywords are found, it appends a consolidated error message to the provided list of errors.
#' This check helps prevent conflicts with internal functionality that relies on these reserved keywords.
#'
#' @param obj The S4 object to check.
#' @param slotName The name of the slot in the object to check for column or vertex names.
#' @param errors A character vector where error messages will be accumulated.
#'               New error messages are appended if reserved keywords are found as column or vertex names in the specified slot.
#'
#' @return A character vector of accumulated error messages, including any new errors found during this check.
#'
#' @importFrom methods slot
#' @importFrom igraph V
#' @keywords internal
.assertNoReservedKeywords <- function(obj, slotName, errors) {
  # Define the list of reserved keywords
  reservedKeywords <- c("sample_id", "taxa_id", "comm_id",
                        "abun", "rela", "norm", 
                        "meta", "taxa", 
                        "netw", "comm",
                        "mgnet", ".")
  
  slotData <- slot(obj, slotName)
  
  # Determine the nature of slotData and set up names to check accordingly
  namesToCheck <- NULL
  if (is.matrix(slotData) || is.data.frame(slotData)) {
    namesToCheck <- colnames(slotData)
  } else if ("igraph" %in% class(slotData)) {
    namesToCheck <- V(slotData)$name
  }
  
  # If applicable, check if names contain reserved keywords
  if (!is.null(namesToCheck)) {
    foundKeywords <- reservedKeywords[reservedKeywords %in% namesToCheck]
    if (length(foundKeywords) > 0) {
      targetType <- ifelse(is.matrix(slotData) || is.data.frame(slotData), "column names", "vertex names")
      message <- sprintf("The following reserved keywords are used as %s in %s slot and cannot be used: {%s}. Please rename these to avoid conflicts with internal functionality.",
                         targetType, slotName, paste(foundKeywords, collapse=", "))
      errors <- c(errors, message)
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
#' @param slotName The name of the slot to compare.
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
  
  netw <- slot(obj, "netw")
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
  if(length(methods::slot(obj, "netw")) == 0){
    errors <- c(errors, "community cannot exist without the associated network.")
  } else if (length(slot(obj, "comm")$membership) != vcount(slot(obj, "netw"))){
    errors <- c(errors, "network and community slots must have the same number of vertices.")
  }
  
  return(errors)
}

#------------------------------------------------------------------------------#
#' Internal: Assert Matching Zero Positions Between Two Slots
#' 
#' Validates that the zero positions in two specified slots within an S4 object are identical.
#' This function is crucial for ensuring data integrity, particularly in biological or ecological
#' datasets where the absence (zero) of a measurement in one type of data must correspond to the absence
#' in another, ensuring aligned analytical outputs.
#'
#' @param obj The S4 object containing the slots to be compared.
#' @param slotName1 The name of the first slot containing a matrix to check for zero positions.
#' @param slotName2 The name of the second slot containing a matrix to check for zero positions.
#' @param errors A pre-existing character vector where any new error messages identified
#'        during the validation process will be appended. This allows for the accumulation
#'        of error messages across multiple validation checks.
#'
#' @return An updated character vector of error messages. If zero positions do not match,
#'         an appropriate error message is added to the vector. If there are no discrepancies,
#'         the original vector is returned unchanged.
#'
#' @importFrom methods slot
#' @keywords internal
.assertMatchingZeroPositions <- function(obj, slotName1, slotName2, errors) {
  
  matrix1 <- slot(obj, slotName1)
  matrix2 <- slot(obj, slotName2)
  
  # Identify mismatching zero positions
  zeroMismatch <- which((matrix1 == 0) != (matrix2 == 0), arr.ind = TRUE)
  
  if (length(zeroMismatch) > 0) {
    # Format the indices into a human-readable string
    errorMessage <- sprintf("Mismatch in zero positions between slots '%s' and '%s.",
                            slotName1, slotName2)
    errors <- c(errors, errorMessage)
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
  errors$abun <- character()
  if( length(object@abun)!=0 ){
    errors$abun <- .assertNumericMatrix(object@abun, errors$abun)
    
    if ( length(errors$abun)==0 ){
      errors$abun <- .assertUniqueRowColNames(object@abun, errors$abun)
      errors$abun <- .assertAllPositive(object@abun, errors$abun)
      errors$abun <- .assertNoReservedKeywords(object, "abun", errors$abun)
    }
  }
  
  #CHECK RELATIVE
  #-------------------------------------#
  errors$rela <- character()
  if( length(object@rela)!=0 ){
    errors$rela <- .assertNumericMatrix(object@rela, errors$rela)
    
    if ( length(errors$rel_rela)==0 ){
      errors$rela <- .assertUniqueRowColNames(object@rela, errors$rela)
      errors$rela <- .assertAllPositive(object@rela, errors$rela)
      errors$rela <- .assertNoReservedKeywords(object, "rela",errors$rela)
    }
  }
  
  #CHECK NORM-ABUNDANCE
  #-------------------------------------#
  errors$norm <- character()
  if( length(object@norm)!=0 ){
    errors$norm <- .assertNumericMatrix(object@norm, errors$norm)
    
    if ( length(errors$norm)==0 ){
      errors$norm <- .assertUniqueRowColNames(object@norm, errors$norm)
      errors$norm <- .assertNoReservedKeywords(object, "norm", errors$norm)
      
    }
  }
  
  #CHECK INFO_SAMPLE
  #-------------------------------------#
  errors$meta <- character()
  if( length(object@meta)!=0 ){
    errors$meta <- .assertDataFrame(object@meta, errors$meta)
    
    if ( length(errors$meta)==0 ){
      errors$meta <- .assertUniqueRowColNames(object@meta, errors$meta)
      errors$meta <- .assertNoReservedKeywords(object, "meta", errors$meta)
    }
  }
  
  #CHECK INFO_TAXA
  #-------------------------------------#
  errors$taxa <- character()
  if( length(object@taxa)!=0 ){
    errors$taxa <- .assertDataFrame(object@taxa, errors$taxa)
    
    if ( length(errors$taxa)==0 ){
      errors$taxa <- .assertUniqueRowColNames(object@taxa, errors$taxa)
      errors$taxa <- .assertNoReservedKeywords(object, "taxa", errors$taxa)
    }
  }
  
  # CHECK NETWORK
  #-------------------------------------#
  errors$netw <- character()
  if(length(object@netw)!=0){
    errors$netw <- .assertNamedIgraph(object@netw, errors$netw)
    errors$netw <- .assertNoReservedKeywords(object, "netw", errors$netw)
  }
  
  # CHECK community
  #-------------------------------------#
  errors$comm <- character()
  if(length(object@comm)!=0){
    errors$comm <- .assertCommunitiesClass(object@comm, errors$comm)
  }
  
  # CHECK RECIPROCAL PROPERTIES
  #-------------------------------------#
  errors$reciprocal <- character()
  
  if(length(object@abun)!=0 && length(errors$abun)==0){
    if(length(object@rela)!=0 && length(errors$rela)==0){
      errors_tmp <- character()
      errors_tmp <- .assertMatchingNames(object, "abun", "rela", errors$errors_tmp)
      if(length(errors_tmp)==0){
        errors$reciprocal <- .assertMatchingZeroPositions(object, "abun", "rela", errors_tmp)}}
    if(length(object@norm)!=0 && length(errors$norm)==0){
      errors$reciprocal <- .assertMatchingNames(object, "abun", "norm", errors$reciprocal)}
    if(length(object@meta)!=0 && length(errors$meta)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "abun", "meta", errors$reciprocal)}
    if(length(object@taxa)!=0 && length(errors$taxa)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "abun", "taxa", errors$reciprocal)}
    if(length(object@netw)!=0 && length(errors$netw)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "abun", "columns", errors$reciprocal)}
  }
  
  if(length(object@rela)!=0 && length(errors$rela)==0){
    if(length(object@norm)!=0 && length(errors$norm)==0){
      errors$reciprocal <- .assertMatchingNames(object, "rela", "norm", errors$reciprocal)}
    if(length(object@meta)!=0 && length(errors$meta)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "rela", "meta", errors$reciprocal)}
    if(length(object@taxa)!=0 && length(errors$taxa)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "rela", "taxa", errors$reciprocal)}
    if(length(object@netw)!=0 && length(errors$netw)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "rela", "columns", errors$reciprocal)}
  }
  
  if(length(object@norm)!=0 && length(errors$norm)==0){
    if(length(object@meta)!=0 && length(errors$meta)==0){
      errors$reciprocal <- .assertMatchingRowNames(object, "norm", "meta", errors$reciprocal)}
    if(length(object@taxa)!=0 && length(errors$taxa)==0){
      errors$reciprocal <- .assertMatchingColsRowsNames(object, "norm", "taxa", errors$reciprocal)}
    if(length(object@netw)!=0 && length(errors$netw)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "norm", "columns", errors$reciprocal)}
  }
  
  if(length(object@taxa)!=0 && length(errors$taxa)==0){
    if(length(object@netw)!=0 && length(errors$netw)==0){
      errors$reciprocal <- .assertMatchingNamesVertices(object, "taxa", "rows", errors$reciprocal)}
  }
  
  if(length(object@comm)!=0){
    errors$reciprocal <- .assertMatchingCommunitiesNetwork(object, errors$reciprocal)
  }
  
  # Consolidate column names from slots they are unique!!
  all_info_names <- character()
  
  if (length(object@taxa) > 0) all_info_names <- c(all_info_names, colnames(object@taxa))
  if (length(object@meta)) all_info_names <- c(all_info_names, colnames(object@meta))
  
  # Add error for duplicated names, if any
  duplicated_names <- unique(all_info_names[which(duplicated(all_info_names))])
  if(length(duplicated_names) > 0) {
    dup_names_str <- paste(duplicated_names, collapse = ", ")
    errors$reciprocal <- c(errors$reciprocal, sprintf("Duplicated column names found across slots: %s. Each column name must be unique across taxa and samples.", dup_names_str))
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