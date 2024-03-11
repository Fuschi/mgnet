#' Internal: Assert Unique Row and Column Names
#'
#' Validates that a given object (matrix or data.frame) has unique row and
#' column names. This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class.
#'
#' @param obj An object to check, either a matrix or a data.frame.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertUniqueRowColNames <- function(obj, slotName) {
  
  errors <- character()
  
  if (!is.null(rownames(obj)) && anyDuplicated(rownames(obj))) {
    errors <- c(errors, sprintf("%s: Row names must be unique.", slotName))
  } else if (is.null(rownames(obj))) {
    errors <- c(errors, sprintf("%s: Object must have row names.", slotName))
  }
  
  if (!is.null(colnames(obj)) && anyDuplicated(colnames(obj))) {
    errors <- c(errors, sprintf("%s: Column names must be unique.", slotName))
  } else if (is.null(colnames(obj))) {
    errors <- c(errors, sprintf("%s: Object must have column names.", slotName))
  }
  return(errors)
}


#' Internal: Assert Abundance Slot
#'
#' Validates that a given matrix is numeric and it has all values >=0. 
#' This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class for the abundance slot.
#'
#' @param obj A matrix to check.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertAbundance <- function(obj, slotName) {
  
  errors <- character()
  errors <- c(errors,.assertUniqueRowColNames(obj))
  
  if (!is.matrix(obj) || !is.numeric(obj)) {
    errors <- c(errors, sprintf("%s must be a numeric matrix", slotName))
  }
  if (any(obj<0, na.rm=TRUE)) {
    errors <- c(errors, sprintf("%s elements must be greater or equal to zero", slotName))
  }
  if (any(is.na(obj))) {
    errors <- c(errors, sprintf("%s cannot be composed from NA elements", slotName))
  }
  
  if(length(errors)>0){
    errors <- c(paste("DEBUGGER of",slotName,"slot"),
                paste("\t-",errors))
  }
  
  return(errors)
}


#' Internal: Assert sampleInfo or taxaInfo Slot
#'
#' Validates that a given data.frame has unique row and column names. 
#' This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class for the sampleInfo slot.
#'
#' @param obj A data.frame to check.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertMetaInfo <- function(obj, slotName) {
  
  errors <- character()
  errors <- c(errors,.assertUniqueRowColNames(obj))
  
  if (!is.data.frame(obj)) {
    errors <- c(errors, sprintf("%s must be a data.frame", slotName))
  }
  
  if(length(errors)>0){
    errors <- c(paste("DEBUGGER of",slotName,"slot"),
                paste("\t-",errors))
  }
  
  return(errors)
}


#' Internal: Assert Lineage Slot
#'
#' Validates that a given matrix is character with unique row and column names. 
#' This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class for the lineage slot.
#'
#' @param obj A matrix to check.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertLineage <- function(obj, slotName) {
  
  errors <- character()
  errors <- c(errors,.assertUniqueRowColNames(obj))
  
  if (!is.matrix(obj) || !is.character(obj)) {
    errors <- c(errors, sprintf("%s must be a character matrix", slotName))
  }
  if (any(is.na(obj))) {
    errors <- c(errors, sprintf("%s cannot be composed from NA elements", slotName))
  }
  
  if(length(errors)>0){
    errors <- c(paste("DEBUGGER of",slotName,"slot"),
                paste("\t-",errors))
  }
  
  return(errors)
}


#' Internal: Assert Association Slot
#'
#' Validates that a given matrix is numeric with unique row and column names. 
#' This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class for the lineage slot.
#'
#' @param obj A matrix to check.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertLineage <- function(obj, slotName) {
  
  errors <- character()
  errors <- c(errors,.assertUniqueRowColNames(obj))
  
  if (!is.matrix(obj) || !is.character(obj)) {
    errors <- c(errors, sprintf("%s must be a character matrix", slotName))
  }
  if (any(is.na(obj))) {
    errors <- c(errors, sprintf("%s cannot be composed from NA elements", slotName))
  }
  
  if(length(errors)>0){
    errors <- c(paste("DEBUGGER of",slotName,"slot"),
                paste("\t-",errors))
  }
  
  return(errors)
}


#' Internal: Assert Association Slot
#'
#' Validates that a given matrix is numeric and square. 
#' This function is designed for internal use within the package
#' to enforce data integrity in the 'mgnet' class for the association slot.
#'
#' @param obj A matrix to check.
#' @param slotName A character string providing the name of the slot being checked,
#' used in the construction of error messages for clarity.
#' @return A character vector containing error messages if the conditions are not met,
#' otherwise an empty character vector indicating success.
#' @keywords internal
#' @noRd
.assertAssociation <- function(obj, slotName) {
  
  errors <- character()
  errors <- c(errors,.assertUniqueRowColNames(obj))
  
  if (!is.matrix(obj) || !is.numeric(obj) || (nrow(obj)!=ncol(obj))) {
    errors <- c(errors, sprintf("%s must be a numeric square matrix", slotName))
  }
  if (any(is.na(obj))) {
    errors <- c(errors, sprintf("%s cannot be composed from NA elements", slotName))
  }
  
  if(length(errors)>0){
    errors <- c(paste("DEBUGGER of",slotName,"slot"),
                paste("\t-",errors))
  }
  
  return(errors)
}