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
#' Return the numeric matrix associated with the ngs data from mg class.
#'
#' @usage data(object)
#'
#' @param object Required \code{\link{mg-class}}.
#'
#' @seealso \code{\link{ngs_data}}
#' @rdname data-methods
#' @docType methods
#' @export
#' @aliases data data
setGeneric("data", function(object) standardGeneric("data"))
#' @rdname data-methods
#' @aliases data,ANY-method
setMethod("data", "ANY", function(object){NULL})
#' @rdname data-methods
#' @aliases data,mg-method
setMethod("data", "mg", function(object){return(object@data)})
#####################################
# META
#####################################
#' Retrieves sample metadata.
#'
#' @description
#' Return the data.frame associated with the sample metadata with experimental
#' variables from mg class.
#'
#' @usage meta(object)
#'
#' @param object Required \code{\link{mg-class}}.
#'
#' @seealso \code{\link{sample_metadata}}
#' @rdname meta-methods
#' @docType methods
#' @export
#' @aliases meta meta
setGeneric("meta", function(object) standardGeneric("meta"))
#' @rdname meta-methods
#' @aliases meta,ANY-method
setMethod("meta", "ANY", function(object){NULL})
#' @rdname meta-methods
#' @aliases meta,mg-method
setMethod("meta", "mg", function(object){return(object@meta)})
#####################################
# TAXA
#####################################
#' Retrieves taxonomy table.
#'
#' @description
#' Return the matrix associated with taxonomy classification from mg class.
#'
#' @usage taxa(object)
#'
#' @param object Required \code{\link{mg-class}}.
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
#' @aliases taxa,mg-method
setMethod("taxa", "mg", function(object){return(object@taxa)})
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
#####################################
#' Assign a new ngs matrix to \code{object}
#'
#' @usage data(object) <- value
#'
#' @param object Required \code{\link{mg-class}}.
#' @param value Required \code{\link{ngs_data-class}} or \code{\link{matrix}}
#'
#' @export
#' @docType methods
#' @rdname assign-ngs_data
#' @aliases assign-ngs_data
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
#' @rdname assign-ngs_data
#' @aliases data<-,mg,ngs_data-method,ngs_data
setMethod("data<-", c("mg", "ngs_data"), function(object, value){
  mg(data=value, meta=object@meta, taxa=object@taxa)
})
#' @rdname assign-ngs_data
#' @aliases data<-,mg,ngs_data-method,matrix
setMethod("data<-", c("mg", "matrix"), function(object, value){
  mg(data=ngs_data(value), meta=object@meta, taxa=object@taxa)
})
#####################################
# META
#####################################
#' Assign new samples metadata to \code{object}
#'
#' @usage meta(object) <- value
#'
#' @param object Required \code{\link{mg-class}}.
#' @param value Required \code{\link{sample_metadata-class}} or \code{\link{data.frame}}
#'
#' @export
#' @docType methods
#' @rdname assign-sample_metadata
#' @aliases assign-sample_metadata
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))
#' @rdname assign-sample_metadata
#' @aliases meta<-,mg,sample_metadata-method,sample_metadata
setMethod("meta<-", c("mg", "sample_metadata"), function(object, value){
  mg(data=object@data, meta=value, taxa=object@taxa)
})
#' @rdname assign-sample_metadata
#' @aliases meta<-,mg,sample_metadata-method,data.frame
setMethod("meta<-", c("mg", "data.frame"), function(object, value){
  mg(data=object@data, meta=sample_metadata(value), taxa=object@taxa)
})
#####################################
# TAXA
#####################################
#' Assign new taxonomy table to \code{object}
#'
#' @usage taxa(object) <- value
#'
#' @param object Required \code{\link{mg-class}}.
#' @param value Required \code{\link{taxonomy_table-class}} or \code{\link{data.frame}}
#'
#' @export
#' @docType methods
#' @rdname assign-taxonomy_table
#' @aliases assign-taxonomy_table
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))
#' @rdname assign-taxonomy_table
#' @aliases taxa<-,mg,taxonomy_table-method,taxonomy_table
setMethod("taxa<-", c("mg", "taxonomy_table"), function(object, value){
  mg(data=object@data, meta=object@meta, taxa=value)
})
#' @rdname assign-taxonomy_table
#' @aliases taxa<-,mg,taxonomy_table-method,matrix
setMethod("taxa<-", c("mg", "matrix"), function(object, value){
  mg(data=object@data, meta=object@meta, taxa=taxonomy_table(value))
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
# DATA
#####################################
#' Method extensions to extraction operator.
#'
#' See the documentation for the \code{\link[base]{Extract}} generic,
#' defined in the R \code{\link[base]{base-package}} for the expected behavior. 
#'
#' @param j See \code{\link[base]{Extract}}
#' 
#' @param ... See \code{\link[base]{Extract}}
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
#' @inheritParams base::Extract
setMethod("[", "ngs_data", function(x, i, j, ...){
  ngs_data(x@.Data[i, j, drop=FALSE, ...])
  })
# META
#####################################
#' @export
#' @rdname extract-methods
setMethod("[", "sample_metadata", function(x, i, j, ...){
  sample_metadata(x@.Data[i, j, drop=FALSE, ...] )
  })
# TAXA
#####################################

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
            cat("*** Class ngs_data , method Show (limited to a matrix 10x10) ***\n")
            nrowShow <- min(10,nrow(object))
            ncolShow <- min(10,ncol(object))
            print(formatC(object[1:nrowShow,1:ncolShow]),quote=FALSE)
            cat("******* End Show (ngs_data) ******* \n")
          })
################################################################################
setMethod("show","sample_metadata",
          function(object){
            cat("*** Class sample_metadata , method Show (limited to a data.frame of 10 rows/samples) ***\n")
            nrowShow <- min(10,nrow(object))
            print(object[1:nrowShow,],quote=FALSE)
            cat("******* End Show (sample_metadata) ******* \n")
          })
################################################################################
setMethod("show","taxonomy_table",
          function(object){
            cat("*** Class taxonomy_table , method Show (limited to a matrix 10x10) *** \n")
            nrowShow <- min(10,nrow(object))
            ncolShow <- min(10,ncol(object))
            print(object[1:nrowShow,1:ncolShow],quote=FALSE)
            cat("******* End Show (taxonomy_table) ******* \n")
          })
################################################################################
setMethod("show","mg",
          function(object){
            cat("*** Class mg , method Show *** \n")
            print(paste("Sample Number:",nrow(data(object))))
            print(paste("Taxa Number:",ncol(data(object))))

            if(is.null(object@meta)){
              print("Sample Meta Data: NA")
            } else {
              print(paste("Sample Meta Data:",paste(colnames(object@meta),collapse="," )))
            }

            if(is.null(object@taxa)){
              print("Taxonomic Ranks: NA")
            } else {
              print(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")))
            }

            cat("******* End Show (mg) ******* \n")
          })
################################################################################
################################################################################
# END SHOW METHODS
################################################################################
################################################################################