################################################################################
################################################################################
# CLASSES
################################################################################
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
setClass(Class="ngs_data",
         contains="matrix", prototype=matrix(nrow=0,ncol=0),
         validity=function(object){
           
           if(any(dim(object)==0)) return("\n matrix must have non-zero dimensions.")
           if(!is.numeric(object)) return("\n matrix must be numeric")
           if(!all(object>=0))return("\n all matrix elements must be greater or equal to zero")
           if(is.null(rownames(object))) return("\n matrix must have the rows names where the samples IDs were saved.")
           if(is.null(colnames(object))) return("\n matrix must have the cols names where the taxa IDs were saved")
           if(any(duplicated(rownames(object)))) return("\n find in matrix at least a duplicated row name / sample ID.")
           if(any(duplicated(colnames(object)))) return("\n find in matrix at least a duplicated col name / taxa ID.")
           
           TRUE
         })
################################################################################
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
setClass("sample_metadata",
         contains="data.frame", prototype=data.frame(),
         validity=function(object){
           
           if(any(dim(object)==0)) return("\n data.frame must have non-zero dimensions.")
           if(is.null(rownames(object))) return("\n matrix must have the rows names where the samples IDs were saved.")
           if(is.null(colnames(object))) return("\n matrix must have the cols names where the experimental variables were saved")
           if(any(duplicated(rownames(object)))) return("\n find in matrix at least a duplicated row name / sample ID.")
           if(any(duplicated(colnames(object)))) return("\n find in matrix at least a duplicated col name / experimental variable.")
           
           TRUE
         })
################################################################################
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
setClass("taxonomy_table",
         contains="matrix", prototype=matrix(nrow=0,ncol=0),
         validity=function(object){
           
           if(any(dim(object)==0)) return("\n matrix must have non-zero dimensions.")
           if(ncol(object)<2) return("\n matrix must have more than a single column/rank. Where are the taxonomic info?")
           if(!is.character(object)) return("\n matrix must be character")
           if(!all(validUTF8(object))) return("\n all matrix elements must be encoded with UTF-8")
           if(is.null(rownames(object))) return("\n matrix must have the rows names where the taxa IDs were saved.")
           if(is.null(colnames(object))) return("\n matrix must have the cols names where the taxonomic ranks were saved")
           if(any(duplicated(rownames(object)))) return("\n find in matrix at least a duplicated row name / taxa ID.")
           if(any(duplicated(colnames(object)))) return("\n find in matrix at least a duplicated col name / rank.")
           if(any(duplicated(rownames(object[,ncol(object)])))) return("\n find in last column matrix at least a duplicated taxa ID.")
           
           # check if the taxonomic hierarchy is well structured. 
           # Basically I don't want two identical taxa that differ from a higher rank.
           # Example:
           # species   genera  family  ...
           #   sp1       g1      f1
           #   sp2       g1      f2
           taxa.rank <- object[,-ncol(object)]
           unique.rank <- unique(taxa.rank[,ncol(taxa.rank)])
           unique.path.rank <- unique(apply(taxa.rank,1,function(x)paste(x,collapse="/")))
           
           if(length(unique.rank)!=length(unique.path.rank))return(paste(
             "\n find error in hierarchy of taxa classification.\n",
             " Two identical taxa can't have different higher taxonomy classification, as example:\n",
             " species   genera  family  ...\n",
             "   sp1       g1      f1    ...\n",
             "   sp2       g1      f2    ...\n",
             " If species sp1 and sp2 belong to genera g1 both can't be classified at family rank with families f1 and f2 differently."
             ,sep=""))
           
           TRUE
         })
################################################################################
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
         taxa="taxonomy_tableOrNULL"),
  
  prototype=prototype(data=NULL,meta=NULL,taxa=NULL),
  
  validity=function(object){
    
    if(!is.null(object@data) && !is.null(object@meta)){
      if(!all(rownames(object@data)==rownames(object@meta))) return("\n rows names / sample IDs must be identical in data and meta slots")
    }
    
    if(!is.null(object@data) && !is.null(object@taxa)){
      if(!all(colnames(object@data)==rownames(object@taxa))) return("\n data colnames and taxa rownames (taxa IDs) must be identical in data and taxa slots")
    }
    
    TRUE
  })
################################################################################
################################################################################
# END CLASSES
################################################################################
################################################################################




################################################################################
################################################################################
# CONSTRUCTORS
################################################################################
################################################################################
#' User constructor for ngs_data s4 class.
#' 
#' User constructor to create an object belonging to formal s4 class ngs_data 
#' avoiding the new function.
#' 
#' @param data numeric matrix with all elements >=0.  
#'@export
ngs_data <- function(data){new("ngs_data", data)}
################################################################################
#' User constructor for sample_metadata s4 class.
#' 
#' User constructor to create an object belonging to formal s4 class sample_metadata 
#' avoiding the new function.
#' 
#' @param meta data.frame with experimental variables.  
#' @export
sample_metadata <- function(meta){new("sample_metadata", meta)}
################################################################################
#' User constructor for taxonomy_table s4 class.
#' 
#' User constructor to create an object belonging to formal s4 class taxonomy_table 
#' avoiding the new function.
#' 
#' @param taxa character matrix with taxonomic classification.  
#' @export
taxonomy_table <- function(taxa){new("taxonomy_table", taxa)}
################################################################################
#' User constructor for mg s4 class.
#' 
#' User constructor to create an object belonging to formal s4 class mg 
#' avoiding the new function.
#' 
#' @param data numeric matrix or ngs_data s4 obj with all elements >=0.  
#' @param meta data.frame or sample_metadata s4 obj with experimental variables.
#' @param taxa character matrix or taxonomy_table s4 obj with taxonomic classification.  
#' @export
mg <- function(data=NULL,meta=NULL,taxa=NULL){
  
  if(class(data)[1]!="ngs_data" && !is.null(data)) data=ngs_data(data)
  if(class(meta)[1]!="sample_metadata" && !is.null(meta)) meta=sample_metadata(meta)
  if(class(taxa)[1]!="taxonomy_table" && !is.null(taxa)) taxa=taxonomy_table(taxa)
  
  return(new("mg",data=data,meta=meta,taxa=taxa))
}
################################################################################
################################################################################
# END CONSTRUCTORS
################################################################################
################################################################################