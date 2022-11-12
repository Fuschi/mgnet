################################################################################
################################################################################
# GETTERS
################################################################################
################################################################################
# DATA
#####################################
#' Retrieves ngs data.
#' 
#' @description 
#' Return the numeric matrix associated with the ngs data from mg or ngs_data 
#' classes.
#'
#' @usage data(object)
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{ngs_data-class}}.
#'
#' @rdname data-methods
#' @docType methods
#' @export
#' @aliases data data
setGeneric("data", function(object) standardGeneric("data"))
#' @rdname data-methods
#' @aliases data,ANY-method
setMethod("data", "ANY", function(object){NULL})
#' @rdname data-methods
#' @aliases data,ngs_data-method
setMethod("data", "ngs_data", function(object){return(object@value)})
#' @rdname data-methods
#' @aliases data,mg-method
setMethod("data", "mg", function(object){
  ifelse(is.null(object@data),return(object@data),return(object@data@value))
})
#####################################
# META
#####################################
#' Retrieves sample metadata.
#'
#' @description
#' Return the data.frame associated with the sample metadata with experimental
#' variables from mg or sample_metadata classes.
#'
#' @usage meta(object)
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{sample_metadata-class}}.
#'
#' @rdname meta-methods
#' @docType methods
#' @export
#' @aliases meta meta
setGeneric("meta", function(object) standardGeneric("meta"))
#' @rdname meta-methods
#' @aliases meta,ANY-method
setMethod("meta", "ANY", function(object){NULL})
#' @rdname meta-methods
#' @aliases meta,sample_metadata-method
setMethod("meta", "sample_metadata", function(object){return(object@value)})
#' @rdname meta-methods
#' @aliases meta,mg-method
setMethod("meta", "mg", function(object){
  ifelse(is.null(object@meta),return(object@meta),return(object@meta@value))
})
#####################################
# TAXA
#####################################
#' Retrieves taxonomy table.
#'
#' @description
#' Return the matrix associated with taxonomy classification from mg or 
#' taxonomy_table classes.
#'
#' @usage taxa(object)
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{taxonomy_table-class}}.
#'
#' @seealso \code{\link{taxonomy_table}}
#' @rdname taxa-methods
#' @docType methods
#' @export
#' @aliases taxa taxa
setGeneric("taxa", function(object) standardGeneric("taxa"))
#' @rdname taxa-methods
#' @aliases taxa,ANY-method
setMethod("taxa", "ANY", function(object){NULL})
#' @rdname taxa-methods
#' @aliases taxa,taxonomy_table-method
setMethod("taxa", "taxonomy_table", function(object){return(object@value)})
#' @rdname taxa-methods
#' @aliases taxa,mg-method
setMethod("taxa", "mg", function(object){
  ifelse(is.null(object@taxa),return(object@taxa),return(object@taxa@value))
})
################################################################################
################################################################################
# END GETTERS
################################################################################
################################################################################




################################################################################
################################################################################
# SETTERS
################################################################################
################################################################################
# DATA
#' @import methods
#' @noRd
setClassUnion("ngs_dataOrMatrixOrNULL", c("ngs_data","matrix", "NULL"))
#####################################
#' Assign a new ngs matrix to \code{object}
#'
#' @usage data(object) <- value
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{ngs_data-class}}.
#' @param value Required \code{\link{ngs_data-class}} or \code{\link{matrix}}
#'
#' @export
#' @docType methods
#' @rdname assign-ngs_data
#' @aliases assign-ngs_data
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
#' @rdname assign-ngs_data
#' @aliases data<-,ngs_data,ngs_data-method
setMethod("data<-", c("ngs_data", "ngs_dataOrMatrixOrNULL"), function(object, value){
  ifelse(class(object)[1]!="ngs_data",return(ngs_data(value)),return(value))
  })
#' @rdname assign-ngs_data
#' @aliases data<-,mg,ngs_data-method
setMethod("data<-", c("mg", "ngs_dataOrMatrixOrNULL"), function(object, value){
  mg(data=value, meta=object@meta, taxa=object@taxa)
})
#####################################
# META
#' @import methods
#' @noRd
setClassUnion("sample_metadataOrDataFrameOrNULL", c("sample_metadata","data.frame", "NULL"))
#####################################
#' Assign new samples metadata to \code{object}
#'
#' @usage meta(object) <- value
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{sample_metadata-class}}.
#' @param value Required \code{\link{sample_metadata-class}} or \code{\link{data.frame}}
#'
#' @export
#' @docType methods
#' @rdname assign-sample_metadata
#' @aliases assign-sample_metadata
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))
#' @rdname assign-sample_metadata
#' @aliases meta<-,sample_metadata,sample_metadata-method
setMethod("meta<-", c("sample_metadata", "sample_metadataOrDataFrameOrNULL"), function(object, value){
  ifelse(class(object)[1]!="sample_metadata",return(sample_metadata(value)),return(value))
})
#' @rdname assign-sample_metadata
#' @aliases meta<-,mg,sample_metadata-method
setMethod("meta<-", c("mg", "ngs_dataOrMatrixOrNULL"), function(object, value){
  mg(data=object@data, meta=value, taxa=object@taxa)
})
#####################################
# TAXA
#' @import methods
#' @noRd
setClassUnion("taxonomy_tableOrMatrixOrNULL", c("taxonomy_table","matrix", "NULL"))
#####################################
#' Assign new taxonomy table to \code{object}
#'
#' @usage taxa(object) <- value
#'
#' @param object Required \code{\link{mg-class}} or \code{\link{taxonomy_table-class}}.
#' @param value Required \code{\link{taxonomy_table-class}} or \code{\link{matrix}}
#'
#' @export
#' @docType methods
#' @rdname assign-taxonomy_table
#' @aliases assign-taxonomy_table
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))
#' @rdname assign-taxonomy_table
#' @aliases taxa<-,taxonomy_table,taxonomy_table-method
setMethod("taxa<-", c("taxonomy_table", "taxonomy_tableOrMatrixOrNULL"), function(object, value){
  ifelse(class(object)[1]!="taxonomy_table",return(taxonomy_table(value)),return(value))
})
#' @rdname assign-taxonomy_table
#' @aliases taxa<-,mg,taxonomy_table-method
setMethod("taxa<-", c("mg", "taxonomy_tableOrMatrixOrNULL"), function(object, value){
  mg(data=object@data, meta=object@meta, taxa=value)
})
################################################################################
################################################################################
# END SETTERS
################################################################################
################################################################################




################################################################################
################################################################################
# EXTRACTORS
################################################################################
################################################################################
#' Method extensions to extraction operator for mgnet objects.
#'
#' @param x See \code{\link[base]{Extract}}
#' @param i See \code{\link[base]{Extract}}
#' @param j See \code{\link[base]{Extract}}
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
setMethod(f="[",
          signature="ngs_data",
          definition=function(x,i,j){
            return(ngs_data(x@value[i,j,drop=FALSE]))
})
#' @export
#' @rdname extract-methods
setMethod(f="[",
          signature="sample_metadata",
          definition=function(x,i,j){
            return(sample_metadata(x@value[i,j,drop=FALSE]))
          })
#' @export
#' @rdname extract-methods
setMethod(f="[",
          signature="taxonomy_table",
          definition=function(x,i,j){
            return(taxonomy_table(x@value[i,j,drop=FALSE]))
          })
#' @export
#' @rdname extract-methods
setMethod(f="[",
          signature="mg",
          definition=function(x,i,j){
            return(mg(data=x@data[i,j,drop=FALSE],
                      meta=x@meta[i, ,drop=FALSE],
                      taxa=x@taxa[ ,j,drop=FALSE]))
          })
################################################################################
################################################################################
# END EXTRACTORS
################################################################################
################################################################################



################################################################################
################################################################################
# SHOW METHODS
################################################################################
################################################################################
setMethod("show","ngs_data",
          function(object){
            
            if(is.null(object@value)){
              cat("*** Class ngs_data , method Show (limited to a matrix 5x5) ***\n")
              cat("NULL")
              cat("******* End Show (ngs_data) ******* \n")
            } else {
              cat("*** Class ngs_data , method Show (limited to a matrix 5x5) ***\n")
              nrowShow <- min(5,nrow(object@value))
              ncolShow <- min(5,ncol(object@value))
              print(formatC(object@value[1:nrowShow,1:ncolShow]),quote=FALSE)
              cat("******* End Show (ngs_data) ******* \n")
            }

          })
################################################################################
setMethod("show","sample_metadata",
          function(object){
            
            if(is.null(object@value)){
              cat("*** Class sample_metadata , method Show (limited to a data.frame of 5 rows/samples) ***\n")
              cat("NULL")
              cat("******* End Show (sample_metadata) ******* \n")
            } else {
              cat("*** Class sample_metadata , method Show (limited to a data.frame of 5 rows/samples) ***\n")
              nrowShow <- min(5,nrow(object@value))
              print(object@value[1:nrowShow,],quote=FALSE)
              cat("******* End Show (sample_metadata) ******* \n")              
            }

          })
################################################################################
setMethod("show","taxonomy_table",
          function(object){
            
            if(is.null(object@value)){
              cat("*** Class taxonomy_table , method Show (limited to 5 taxa) *** \n")
              cat("NULL")
              cat("******* End Show (taxonomy_table) ******* \n")
            } else {
              cat("*** Class taxonomy_table , method Show (limited to 5 taxa) *** \n")
              nrowShow <- min(5,nrow(object@value))
              print(object@value[1:nrowShow,],quote=FALSE)
              cat("******* End Show (taxonomy_table) ******* \n")
            }

          })
################################################################################
setMethod("show","mg",
          function(object){
            cat("*** Class mg , method Show *** \n")
            cat(paste("Sample Number:",max(nrow(object@data@value),nrow(object@meta@value)),"\n"))
            cat(paste("Taxa Number:",max(ncol(object@data@value),nrow(object@taxa@value)),"\n"))
            
            if(is.null(object@meta)){
              cat("Sample Meta Data: NA \n")
            } else {
              cat(paste("Sample Meta Data:",paste(colnames(object@meta@value),collapse="," )),"\n")
            }
            
            if(is.null(object@taxa)){
              cat("Taxonomic Ranks: NA")
            } else {
              cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa@value),collapse=",")),"\n")
            }
            cat("\n")
            print(object@data)
            cat("\n")
            print(object@meta)
            cat("\n")
            print(object@taxa)
            cat("\n")
            cat("******* End Show (mg) ******* \n")
          })
################################################################################
################################################################################
# END SHOW METHODS
################################################################################
################################################################################