################################################################################
#' The S4 class to store NGS data.
#'
#' @description
#' The information is saved in a  \code{\link{matrix}} where the rows express the samples and 
#' the columns the taxa.
#'
#' \describe{
#'    \item{.Data}{This slot is the \code{\link{matrix}} containing NGS data.}
#'  }
#'
#' @name NGS_data-class
#' @rdname NGS_data-class
#' @exportClass NGS_data
NGS_data <- setClass("NGS_data", contains="matrix")
################################################################################


#' @noRd
setClassUnion("NGS_dataOrNULL", c("NGS_data", "NULL"))
#' @noRd
setClassUnion("MatrixOrNULL", c("matrix", "NULL"))
#' @noRd 
setClassUnion("DataFrameOrNULL", c("data.frame", "NULL"))



#' S4 class to manage metagenomics data
#'
#'
#' @slot data numeric matrix with data
#' @slot meta sample features as data.frame
#' @slot taxa character matrix with taxonomic information 
#'
#' @name mg-class
#' @rdname mg-class
#' @exportClass mg
mg <- setClass(
  Class="mg",
  
  slot=c(data="NGS_data",
         meta="DataFrameOrNULL",
         taxa="MatrixOrNULL")
)