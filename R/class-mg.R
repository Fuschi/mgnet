################################################################################
################################################################################
# CLASS MG
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
  
  slot=c(data="matrix",
         meta="data.frame",
         taxa="matrix"),
  
  prototype=prototype(data=matrix(nrow=0,ncol=0),
                      meta=data.frame(),
                      taxa=matrix(nrow=0,ncol=0)),
  
  validity=function(object){
    
    #CHECK DATA
    #-------------------------------------#
    if( length(object@data)!=0 ){
      if(!is.numeric(object@data)) return("\n data matrix must be numeric")
      if(!all(object@data>=0))return("\n all data matrix elements must be greater or equal to zero")
      if(is.null(rownames(object@data))) return("\n data matrix must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@data))) return("\n data matrix must have the cols names where the taxa IDs were saved")
      if(any(duplicated(rownames(object@data)))) return("\n find in data matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@data)))) return("\n find in data matrix at least a duplicated col name / taxa ID.")
    }
    
    #CHECK META
    #-------------------------------------#
    if( length(object@meta)!=0 ){
      if(is.null(rownames(object@meta))) return("\n meta matrix must have the rows names where the samples IDs were saved.")
      if(is.null(colnames(object@meta))) return("\n meta matrix must have the cols names where the experimental variables were saved")
      if(any(duplicated(rownames(object@meta)))) return("\n find in meta matrix at least a duplicated row name / sample ID.")
      if(any(duplicated(colnames(object@meta)))) return("\n find in meta matrix at least a duplicated col name / experimental variable.")
    }
    
    #CHECK TAXA
    #-------------------------------------#
    if( length(object@taxa)!=0 ){
      if(!is.character(object@taxa)) return("\n taxa matrix must be character")
      if(!all(validUTF8(object@taxa))) return("\n all taxa matrix elements must be encoded with UTF-8")
      if(is.null(rownames(object@taxa))) return("\n taxa matrix must have the rows names where the taxa IDs were saved.")
      if(is.null(colnames(object@taxa))) return("\n taxa matrix must have the cols names where the taxonomic ranks were saved")
      if(any(duplicated(rownames(object@taxa)))) return("\n find in taxa matrix at least a duplicated row name / taxa ID.")
      if(any(duplicated(colnames(object@taxa)))) return("\n find in taxa matrix at least a duplicated col name / rank.")
      if(any(duplicated(rownames(object@taxa[,ncol(object@taxa)])))) return("\n find in last column taxa matrix at least a duplicated taxa ID.")
      
      # check if the taxonomic hierarchy is well structured. 
      if(ncol(object@taxa)>=2){
        taxa.rank <- object@taxa[,-ncol(object@taxa),drop=FALSE]
        unique.rank <- unique(taxa.rank[,ncol(taxa.rank),drop=FALSE])
        unique.path.rank <- unique(apply(taxa.rank,1,function(x)paste(x,collapse="/")))
        
        if(length(unique.rank)!=length(unique.path.rank))return(paste(
          "\n find error in hierarchy of taxa classification.\n",
          " Two identical taxa can't have different higher taxonomy classification, as example:\n",
          " species   genera  family  ...\n",
          "   sp1       g1      f1    ...\n",
          "   sp2       g1      f2    ...\n",
          " If species sp1 and sp2 belong to genera g1 both can't be classified at family rank with families f1 and f2 differently."
          ,sep=""))
      }}
    
    #CHECK SAMPLE/TAXA ORDER IN SOLTS
    #-------------------------------------#
    if(length(object@data)!=0 && length(object@meta)!=0){
      if(nrow(object@data)!=nrow(object@meta)) return("\n different number of samples in data and meta slots")
      if(!all(rownames(object@data)==rownames(object@meta))) return("\n rows names / sample IDs must be identical in data and meta slots")
    }
    
    if(length(object@data)!=0 && length(object@taxa)!=0){
      if(ncol(object@data)!=nrow(object@taxa)) return("\n different number of taxa in data and taxa slots")
      if(!all(colnames(object@data)==rownames(object@taxa))) return("\n data colnames and taxa rownames (taxa IDs) must be identical in data and taxa slots")
    }
    
    TRUE
  })
################################################################################
################################################################################
# END CLASS MG
################################################################################
################################################################################




################################################################################
################################################################################
# CONSTRUCTOR MG
################################################################################
################################################################################
#' User constructor for mg s4 class.
#' 
#' User constructor to create an object belonging to formal s4 class mg 
#' avoiding the new function.
#' 
#' @param data numeric matrix with all elements >=0.  
#' @param meta data.frame with experimental variables.
#' @param taxa character matrix with taxonomic classification.  
#' @export
mg <- function(data=matrix(nrow=0,ncol=0),
               meta=data.frame(),
               taxa=matrix(nrow=0,ncol=0)){
  
  return(new("mg",data=data,meta=meta,taxa=taxa))
}
################################################################################
################################################################################
# END CONSTRUCTOR MG
################################################################################
################################################################################




################################################################################
################################################################################
# GETTERS MG
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
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname data-methods
#' @docType methods
#' @export
#' @aliases data data
setGeneric("data", function(object) standardGeneric("data"))
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
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname meta-methods
#' @docType methods
#' @export
#' @aliases meta meta
setGeneric("meta", function(object) standardGeneric("meta"))
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
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname taxa-methods
#' @docType methods
#' @export
#' @aliases taxa taxa
setGeneric("taxa", function(object) standardGeneric("taxa"))
#' @rdname taxa-methods
#' @aliases taxa,mg-method
setMethod("taxa", "mg", function(object){return(object@taxa)})
################################################################################
################################################################################
# END GETTERS MG
################################################################################
################################################################################




################################################################################
################################################################################
# SETTERS MG
################################################################################
################################################################################
# DATA<-
#####################################
#' Assign a new ngs matrix to \code{object}
#'
#' @usage data(object) <- value
#'
#' @param object (Required) \code{\link{mg-class}}.
#' @param value (Required) \code{\link{matrix}}
#'
#' @export
#' @docType methods
#' @rdname assign-data
#' @aliases assign-data
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
#' @rdname assign-data
#' @aliases data<-,mg,data-method
setMethod("data<-", c("mg", "matrix"), function(object, value){
  new("mg",data=value, meta=object@meta, taxa=object@taxa)
})
#####################################
# META<-
#####################################
#' Assign new samples metadata to \code{object}
#'
#' @usage meta(object) <- value
#'
#' @param object (Required) \code{\link{mg-class}}.
#' @param value (Required) \code{\link{data.frame}}
#'
#' @export
#' @docType methods
#' @rdname assign-meta
#' @aliases assign-meta
setGeneric("meta<-", function(object, value) standardGeneric("meta<-"))
#' @rdname assign-meta
#' @aliases meta<-,mg,meta-method
setMethod("meta<-", c("mg", "data.frame"), function(object, value){
  new("mg",data=object@data, meta=value, taxa=object@taxa)
})
#####################################
# TAXA<-
#####################################
#' Assign new taxonomy table to \code{object}
#'
#' @usage taxa(object) <- value
#'
#' @param object (Required) \code{\link{mg-class}}.
#' @param value (Required) \code{\link{matrix}}
#'
#' @export
#' @docType methods
#' @rdname assign-taxa
#' @aliases assign-taxa
setGeneric("taxa<-", function(object, value) standardGeneric("taxa<-"))
#' @rdname assign-taxa
#' @aliases taxa<-,mg,taxa-method
setMethod("taxa<-", c("mg", "matrix"), function(object, value){
  new("mg",data=object@data, meta=object@meta, taxa=value)
})
################################################################################
################################################################################
# END SETTERS MG
################################################################################
################################################################################



################################################################################
################################################################################
# EXTRACTOR MG
################################################################################
################################################################################
#' Method extensions to extraction operator for mg object.
#'
#' @param x See \code{\link[base]{Extract}}, \code{\link{mg}} object.
#' @param i See \code{\link[base]{Extract}}, samples indices.
#' @param j See \code{\link[base]{Extract}}, taxa indices
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' 
#' @rdname extract-methods
setMethod(f="[",
          signature="mg",
          definition=function(x,i,j){
            return(new("mg",
                data=x@data[i,j,drop=FALSE],
                meta=x@meta[i, ,drop=FALSE],
                taxa=x@taxa[j, ,drop=FALSE]))
          })
################################################################################
################################################################################
# END EXTRACTOR MG
################################################################################
################################################################################




################################################################################
################################################################################
# SHOW METHOD MG
################################################################################
################################################################################
setMethod("show","mg",
          function(object){
            cat("******* Class mg , method Show ******* \n")
            cat(paste("Sample Number:",max(nrow(object@data),nrow(object@meta)),"\n"))
            cat(paste("Taxa Number:",max(ncol(object@data),nrow(object@taxa)),"\n"))
            cat(paste("Sample Meta Data:",paste(colnames(object@meta),collapse="," )),"\n")
            cat(paste("Taxonomic Ranks:",paste(colnames(object@taxa),collapse=",")),"\n")
            cat("@data (limited to a 3x3 matrix) \n")
            if(length(object@data!=0)){
              nrowShow <- min(3,nrow(object@data))
              ncolShow <- min(3,ncol(object@data))
              print(object@data[1:nrowShow,1:ncolShow],quote=FALSE)
            } else {
              print(object@data)
            }
            cat("@meta (limited to a 3x3 data.frame) \n")
            if(length(object@meta!=0)){
              nrowShow <- min(3,nrow(object@meta))
              ncolShow <- min(3,ncol(object@meta))
              print(object@meta[1:nrowShow,1:ncolShow],quote=FALSE)
            } else {
              print(object@meta)
            }
            cat("@taxa (limited to a 3x3 matrix) \n")
            if(length(object@taxa)!=0){
              nrowShow <- min(3,nrow(object@taxa))
              ncolShow <- min(3,ncol(object@taxa))
              print(object@taxa[1:nrowShow,1:ncolShow],quote=FALSE)
            } else {
              print(object@taxa)
            }
            cat("********** End Show (mg) ********** \n")
          })
################################################################################
################################################################################
# END SHOW METHOD MG
################################################################################
################################################################################




################################################################################
################################################################################
# BASE METHODS
################################################################################
################################################################################
# NSAMPLE
#####################################
#' Get number of samples.
#' 
#' @description 
#' Return an integer indicating the number of sample.
#'
#' @usage nsample(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname nsample-methods
#' @docType methods
#' @export
#' @aliases nsample nsample
setGeneric("nsample", function(object) standardGeneric("nsample"))
#' @rdname nsample-methods
#' @aliases nsample,mg-method
setMethod("nsample", "mg", function(object){
  if(length(object@data!=0)) return(nrow(object@data))
  else if(length(object@meta!=0)) return(nrow(object@meta))
  else return(0)
})
#####################################
# NTAXA
#####################################
#' Get number of taxa
#' 
#' @description 
#' Return an integer indicating the number of taxa
#'
#' @usage ntaxa(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname ntaxa-methods
#' @docType methods
#' @export
#' @aliases ntaxa ntaxa
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))
#' @rdname ntaxa-methods
#' @aliases ntaxa,mg-method
setMethod("ntaxa", "mg", function(object){
  if(length(object@data!=0)) return(ncol(object@data))
  else if(length(object@taxa!=0)) return(nrow(object@taxa))
  else return(0)
})
#####################################
# SAMPLE NAME
#####################################
#' Get samples names.
#' 
#' @description 
#' Return names of samples as character vector.
#'
#' @usage sample_name(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname sample_name-methods
#' @docType methods
#' @export
#' @aliases sample_name sample_name
setGeneric("sample_name", function(object) standardGeneric("sample_name"))
#' @rdname sample_name-methods
#' @aliases sample_name,mg-method
setMethod("sample_name", "mg", function(object){
  if(length(object@data!=0)) return(rownames(object@data))
  else if(length(object@taxa!=0)) return(rownames(object@meta))
  else return(0)
})
#####################################
# TAXA ID
#####################################
#' Get taxa ID.
#' 
#' @description 
#' Return taxonomy id as a character vector.
#'
#' @usage taxaID(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname taxaID-methods
#' @docType methods
#' @export
#' @aliases taxaID taxaID
setGeneric("taxaID", function(object) standardGeneric("taxaID"))
#' @rdname taxaID-methods
#' @aliases taxaID,mg-method
setMethod("taxaID", "mg", function(object){
  if(length(object@data!=0)) return(colnames(object@data))
  else if(length(object@taxa!=0)) return(rownames(object@taxa))
  else return(0)
})
#####################################
# RANKS 
#####################################
#' Get taxonomic ranks.
#' 
#' @description 
#' Return taxonomy ranks as a character vector.
#'
#' @usage ranks(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname ranks-methods
#' @docType methods
#' @export
#' @aliases ranks ranks
setGeneric("ranks", function(object) standardGeneric("ranks"))
#' @rdname ranks-methods
#' @aliases ranks,mg-method
setMethod("ranks", "mg", function(object){return(colnames(object@taxa))})
#####################################
# NRANKS 
#####################################
#' Get taxonomic ranks number.
#' 
#' @description 
#' Return taxonomy ranks number number as integer.
#'
#' @usage nrank(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname nrank-methods
#' @docType methods
#' @export
#' @aliases nrank nrank
setGeneric("nrank", function(object) standardGeneric("nrank"))
#' @rdname nrank-methods
#' @aliases nrank,mg-method
setMethod("nrank", "mg", function(object){return(ncol(object@taxa))})
#####################################
# SAMPLE_INFO 
#####################################
#' Get sample metadata variables.
#' 
#' @description 
#' Return sample metadata variables number as character vector.
#'
#' @usage sample_info(object)
#'
#' @param object (Required) \code{\link{mg-class}}.
#'
#' @rdname sample_info-methods
#' @docType methods
#' @export
#' @aliases sample_info sample_info
setGeneric("sample_info", function(object) standardGeneric("sample_info"))
#' @rdname sample_info-methods
#' @aliases sample_info,mg-method
setMethod("sample_info", "mg", function(object){return(colnames(object@meta))})
#####################################
# TAXA NAME 
#####################################
#' Get taxa name.
#' 
#' @description 
#' Get taxa name at choosen rank as character vector.
#' 
#' @usage taxa_name(object, rank)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param rank taxonomic level choosen (if it is not set, the finest taxonomic rank is assumed)
#' 
#' @rdname taxa_name-methods
#' @docType methods
#' @export
#' @aliases taxa_name taxa_name
setGeneric("taxa_name", function(object, rank) standardGeneric("taxa_name"))
#' @rdname taxa_name-methods
#' @aliases taxa_name,mg-method,missing
setMethod("taxa_name", c("mg","missing"), function(object) object@taxa[,nrank(object)])
#' @rdname taxa_name-methods
#' @aliases taxa_name,mg-method,character
setMethod("taxa_name", c("mg","character"),
          function(object, rank){
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",
                                                    toString(ranks(object)),"}"))
            return(object@taxa[,rank])
          })
################################################################################
################################################################################
# END BASE METHODS
################################################################################
################################################################################




################################################################################
################################################################################
# AGGREGATE TAXA
################################################################################
################################################################################
#' Organize data in higher taxonomic level.
#' 
#' @description 
#' Reorganize an \code{\link{mg-class}} object in an higher taxonomy rank.
#' The function sums the taxa with the same classification and return a new mg
#' object with appropriate data and taxa slots.
#' 
#' @usage aggregate_taxa(object, rank)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param rank taxonomic level choosen.
#' 
#' @rdname aggregate_taxa-methods
#' @docType methods
#' @export
#' @aliases aggregate_taxa aggregate_taxa
#' @export
setGeneric("aggregate_taxa", function(object, rank) standardGeneric("aggregate_taxa"))
#' @rdname taxa_name-methods
#' @aliases aggregate_taxa,mg-method,character
setMethod("aggregate_taxa", c("mg","character"),
          function(object, rank){
            
            if(length(object@data)==0 || length(object@data)==0) stop("data and taxa slots must be present")
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
            
            different.taxa <- unique(object@taxa[,rank])
            data.aggregate <- data.frame(matrix(NA, nrow=nsample(object), ncol=length(different.taxa),
                                                dimnames=list(sample_name(object),different.taxa)))
            
            for(taxa.i in different.taxa){
              idx <- which(taxa.i == object@taxa[,rank])
              data.aggregate[,taxa.i] <- apply(X=object@data, MARGIN=1, function(x) sum(x[idx]) )
            }
            
            taxa.aggregate <- object@taxa[,1:which(ranks(object)==rank)]
            taxa.aggregate <- taxa.aggregate[!duplicated(taxa.aggregate), ]
            rownames(taxa.aggregate) <- taxa.aggregate[,rank]
            taxa.aggregate <- taxa.aggregate[colnames(data.aggregate),]
            
            return(new("mg",data=data.aggregate, meta=object@meta, taxa=taxa.aggregate))
          })
################################################################################
################################################################################
# END AGGREGATE TAXA
################################################################################
################################################################################




################################################################################
################################################################################
# MGMELT
################################################################################
################################################################################
#' Melt mg object into form suitable for easy casting.
#' 
#' @description 
#' Summarize all elements present in an mg object in a single data frame. 
#' The functioning is similar to melt function of reshape2 package and it will 
#' be useful as a preprocess for ggplot2.
#' 
#' @usage mgmelt(object)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' 
#' @importFrom reshape2 melt
#' @rdname mgmelt-methods
#' @docType methods
#' @export
#' @aliases mgmelt mgmelt
#' @export
setGeneric("mgmelt", function(object) standardGeneric("mgmelt"))
#' @rdname mgmelt-methods
#' @aliases mgmelt,mg-method
setMethod("mgmelt", "mg",
          function(object){
            
            if(length(object@data)==0){stop("\n data slot must be present")}
            
            mdf <- melt(data=object@data)
            colnames(mdf) <- c("SampleID","TaxaID","Abundance")
            rownames(mdf) <- paste(mdf$SampleID,"-",mdf$TaxaID,sep="")
            
            if(length(object@taxa!=0)) mdf <- cbind(mdf,object@taxa[mdf$TaxaID,])
            if(length(object@meta!=0)) mdf <- cbind(mdf,object@meta[mdf$SampleID,])
            mdf <- mdf[,-which(duplicated(t(mdf)))]
            
            return(mdf)
          })
################################################################################
################################################################################
# END MGMELT
################################################################################
################################################################################




################################################################################
################################################################################
# AGGREGATE TAXA
################################################################################
################################################################################
#' Organize data in higher taxonomic level.
#' 
#' @description 
#' This provides a convenient way to aggregate mg object taxa. Calculates the
#' sum of abundances over all taxa that map to the same higher-level group.
#' 
#' @usage aggregate_taxa(object,rank)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param rank (Required) character indicates the taxonomic level choosen.
#' 
#' @rdname aggregate_taxa-methods
#' @docType methods
#' @export
#' @aliases aggregate_taxa aggregate_taxa
#' @export
setGeneric("aggregate_taxa", function(object, rank) standardGeneric("aggregate_taxa"))
#' @rdname aggregate_taxa-methods
#' @aliases aggregate_taxa,mg-method
setMethod("aggregate_taxa", c("mg","character"),
          function(object, rank){
            
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
            
            different.taxa <- unique(object@taxa[,rank])
            data.aggregate <- matrix(NA, nrow=nsample(object), ncol=length(different.taxa),
                                                dimnames=list(sample_name(object),different.taxa))
            
            for(taxa.i in different.taxa){
              idx <- which(taxa.i == object@taxa[,rank])
              data.aggregate[,taxa.i] <- apply(X=object@data, MARGIN=1, function(x) sum(x[idx]) )
            }
            
            taxa.aggregate <- object@taxa[,1:which(ranks(object)==rank)]
            taxa.aggregate <- taxa.aggregate[!duplicated(taxa.aggregate), ]
            rownames(taxa.aggregate) <- taxa.aggregate[,rank]
            taxa.aggregate <- taxa.aggregate[colnames(data.aggregate),]
            
            return(new("mg",data=data.aggregate, meta=object@meta, taxa=taxa.aggregate))
          })
################################################################################
################################################################################
# END MGMELT
################################################################################
################################################################################