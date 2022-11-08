################################################################################
#' The S4 class to store ngs data.
#'
#' @description
#' The information is saved in a  \code{\link{matrix}} where the rows express the samples and 
#' the columns the taxa.
#'
#' \describe{
#'    \item{.Data}{Slot containing \code{\link{matrix}} with ngs data.}
#'  }
#'
#' @import methods
#' @name ngs_data-class
#' @rdname ngs_data-class
#' @exportClass ngs_data
setClass("ngs_data", contains="matrix")
################################################################################
#' The S4 class for storing sample metadata
#'
#' @description
#' Rows represent samples, while column indices represent experimental variables
#' that describe the samples.
#'
#' \describe{
#'    \item{.Data}{Slot containing \code{\link{data.frame}} with sample metadata.}
#'  }
#' 
#' @import methods
#' @name sample_metadata-class
#' @rdname sample_metadata-class
#' @exportClass sample_metadata
setClass("sample_metadata", contains="data.frame")
################################################################################
#' An S4 class for storing taxonomic classification of taxa as a character
#' matrix.
#'
#' Row indices represent taxa, columns represent taxonomic ranks.
#' 
#' \describe{
#'    \item{.Data}{Slot containing \code{\link{matrix}} class.}
#' }
#'
#' @import methods
#' @name taxonomy_table-class
#' @rdname taxonomy_table-class
#' @exportClass taxonomy_table
setClass("taxonomy_table", contains="matrix")
################################################################################
#' @import methods
#' @noRd
setClassUnion("ngs_dataOrNULL", c("ngs_data", "NULL"))
#' @import methods
#' @noRd
setClassUnion("sample_metadataOrNULL", c("sample_metadata", "NULL"))
#' @import methods
#' @noRd 
setClassUnion("taxonomy_tableOrNULL", c("taxonomy_table", "NULL"))
################################################################################
#' S4 class to manage metagenomics data
#'
#'
#' @slot data ngs data
#' @slot meta samples experimental features.
#' @slot taxa taxonomic table 
#'
#' @import methods
#' @name mg-class
#' @rdname mg-class
#' @exportClass mg
setClass(
  Class="mg",
  
  slot=c(data="ngs_dataOrNULL",
         meta="sample_metadataOrNULL",
         taxa="taxonomy_tableOrNULL")
)