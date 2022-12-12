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
  
  if(!("Depth" %in% colnames(meta)) && length(data)!=0){
    cat("******* class mg constructor message *******\n")
    cat("Add the samples depth as additional column named Depth to meta\n")
    cat("(Depths are calculated summing counts on each row/sample of data slot)\n")
    cat("It can be accessed with the depth method\n")
    ifelse(length(meta)!=0,
           meta$Depth <- rowSums(data),
           meta <- data.frame("Depth"=rowSums(data)))
    cat("********************************************\n")
  }
  
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
#' @rdname data
#' @docType methods
#' @export
setGeneric("data", function(object) standardGeneric("data"))
#' @rdname data
#' @aliases data,mg
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
#' @rdname meta
#' @docType methods
#' @export
setGeneric("meta", function(object) standardGeneric("meta"))
#' @rdname meta
#' @aliases meta,mg
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
#' @rdname taxa
#' @docType methods
#' @export
setGeneric("taxa", function(object) standardGeneric("taxa"))
#' @rdname taxa
#' @aliases taxa,mg
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
setGeneric("data<-", function(object, value) standardGeneric("data<-"))
#' @rdname assign-data
#' @aliases data<-,mg,matrix
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
#' @aliases meta<-,mg,data.frame
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
#' @aliases taxa<-,mg,matrix
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
#' @rdname nsample
#' @docType methods
#' @export
setGeneric("nsample", function(object) standardGeneric("nsample"))
#' @rdname nsample
#' @aliases nsample,mg
setMethod("nsample", "mg", function(object){
  if(length(object@data!=0)) return(nrow(object@data))
  else if(length(object@meta!=0)) return(nrow(object@meta))
  else return(NULL)
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
#' @rdname ntaxa
#' @docType methods
#' @export
setGeneric("ntaxa", function(object) standardGeneric("ntaxa"))
#' @rdname ntaxa
#' @aliases ntaxa,mg
setMethod("ntaxa", "mg", function(object){
  if(length(object@data!=0)) return(ncol(object@data))
  else if(length(object@taxa!=0)) return(nrow(object@taxa))
  else return(NULL)
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
#' @rdname sample_name
#' @docType methods
#' @export
setGeneric("sample_name", function(object) standardGeneric("sample_name"))
#' @rdname sample_name
#' @aliases sample_name,mg
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
#' @rdname taxaID
#' @docType methods
#' @export
setGeneric("taxaID", function(object) standardGeneric("taxaID"))
#' @rdname taxaID
#' @aliases taxaID,mg
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
#' @rdname ranks
#' @docType methods
#' @export
setGeneric("ranks", function(object) standardGeneric("ranks"))
#' @rdname ranks
#' @aliases ranks,mg
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
#' @rdname nrank
#' @docType methods
#' @export
setGeneric("nrank", function(object) standardGeneric("nrank"))
#' @rdname nrank
#' @aliases nrank,mg
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
#' @rdname sample_info
#' @docType methods
#' @export
setGeneric("sample_info", function(object) standardGeneric("sample_info"))
#' @rdname sample_info
#' @aliases sample_info,mg
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
#' @rdname taxa_name
#' @docType methods
#' @export
setGeneric("taxa_name", function(object, rank) standardGeneric("taxa_name"))
#' @rdname taxa_name
#' @aliases taxa_name,mg,missing
setMethod("taxa_name", c("mg","missing"), function(object) object@taxa[,nrank(object)])
#' @rdname taxa_name
#' @aliases taxa_name,mg,character
setMethod("taxa_name", c("mg","character"),
          function(object, rank){
            if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",
                                                    toString(ranks(object)),"}"))
            return(object@taxa[,rank])
          })
#####################################
# DEPTH 
#####################################
#' Get samples depth.
#' 
#' @description 
#' Retrieves the depth of each sample using the Depth variable in in meta slot.
#' Depth is automatically generated if not present during the construction of 
#' the mg object.
#' 
#' @usage depth(object)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' 
#' @rdname depth
#' @docType methods
#' @export
setGeneric("depth", function(object) standardGeneric("depth"))
#' @rdname depth
#' @aliases depth,mg
setMethod("depth", "mg",function(object)return(object@meta$Depth))
#####################################
# ABUNDANCE 
#####################################
#' Get abundances at choosen rank.
#' 
#' @description 
#' Retrieves the abundances of data at choosen taxonomy rank. Abundance at 
#' higher rank are returned as sums of abundance of elements with the same
#' classification.
#' 
#' @usage abundance(object,rank)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param rank (Optional) character with the taxonomic rank choosen.
#' 
#' @export
#' @docType methods
#' @rdname abundance
setGeneric("abundance", function(object,rank) standardGeneric("abundance"))
#' @rdname abundance
#' @aliases abundance,missing
setMethod("abundance", c("mg","missing"),function(object)return(object@data))
#' @rdname abundance
#' @aliases abundance,character
setMethod("abundance", c("mg","character"),function(object,rank){
  
  if(length(object@data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
  if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
  
  different.taxa <- unique(object@taxa[,rank])
  data.aggregate <- matrix(NA, nrow=nsample(object), ncol=length(different.taxa),
                           dimnames=list(sample_name(object),different.taxa))
  
  for(taxa.i in different.taxa){
    idx <- which(taxa.i == object@taxa[,rank])
    data.aggregate[,taxa.i] <- apply(X=object@data, MARGIN=1, function(x) sum(x[idx]) )
  }
  
  return(data.aggregate)
})
#####################################
# RELVATIVE 
#####################################
#' Get relative abundances.
#' 
#' @description 
#' Retrieves the relative abundances of data, normalized by their depths at
#' the taxonomic rank choosen.
#' 
#' @usage relative(object,rank)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param rank (Optional)
#' 
#' @rdname relative
#' @docType methods
#' @export
setGeneric("relative", function(object,rank) standardGeneric("relative"))
#' @rdname relative
#' @aliases relative,mg,missing
setMethod("relative", c("mg","missing"),function(object)return(object@data/depth(object)))
#' @rdname relative
#' @aliases relative,mg,character
setMethod("relative", c("mg","character"),function(object,rank){
  if(length(object@data)==0 || length(object@taxa)==0) stop("data and taxa slots must be present")
  if(!(rank%in%ranks(object))) stop(paste("rank must be one this possible choises {",toString(ranks(object)),"}"))
  
  different.taxa <- unique(object@taxa[,rank])
  data.aggregate <- matrix(NA, nrow=nsample(object), ncol=length(different.taxa),
                           dimnames=list(sample_name(object),different.taxa))
  
  for(taxa.i in different.taxa){
    idx <- which(taxa.i == object@taxa[,rank])
    data.aggregate[,taxa.i] <- apply(X=object@data, MARGIN=1, function(x) sum(x[idx]) )
  }
  
  return(data.aggregate/depth(object))
})
#####################################
# EMPTY 
#####################################
#' Check if the object is empty.
#' 
#' @description 
#' Control if all slots have lengths equal to 0 (therefore empty).
#' 
#' @usage empty(object)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' 
#' @rdname empty
#' @docType methods
#' @export
setGeneric("empty", function(object) standardGeneric("empty"))
#' @rdname empty
#' @aliases empty,mg
setMethod("empty", c("mg"),function(object){
  length(object@data)==0 & length(object@meta)==0 & length(object@taxa)==0
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
#' @rdname aggregate_taxa
#' @docType methods
#' @export
#' @export
setGeneric("aggregate_taxa", function(object, rank) standardGeneric("aggregate_taxa"))
#' @rdname aggregate_taxa
#' @aliases aggregate_taxa,mg,character
setMethod("aggregate_taxa", c("mg","character"),
          function(object, rank){
            
            if(length(object@data)==0 || length(object@data)==0) stop("data and taxa slots must be present")
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
            mdf$Relative <- melt(relative(object))$value
            
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
# FILTER TAXA
################################################################################
################################################################################
#' Filter taxa based on across-sample abundance criteria.
#' 
#' @description 
#' It applies an arbitrary set of functions list as across-sample criteria,
#' one taxa at a time. The function takes as input a mg object,
#' and returns its trimmed mg version or a logical vector
#' indicating whether or not each taxa passed the criteria.
#' Alternatively, if the \code{"trim"} option is set to \code{FALSE},
#' it return a logical vector indicating whether or not each taxa passed 
#' the criteria.
#' 
#' @usage filter_taxa(object, flist, join.trim)
#' 
#' @param object (Required) \code{\link{mg-class}}.
#' @param flist (Required) \code{\link{list}}. Each element of flist it must be
#' a function.
#' @param join.trim (Optional) Default \code{FALSE}.
#' 
#' @rdname filter_taxa-methods
#' @docType methods
#' @export
setGeneric("filter_taxa", function(object,flist,join.trim) standardGeneric("filter_taxa"))
#' @rdname filter_taxa-methods
#' @aliases filter_taxa,mg,logical
setMethod("filter_taxa", c("mg","list","logical"),
          function(object,flist,join.trim){
            
            if( any(unlist(lapply(flist,class))!="function") ){stop("all flist elements must be a function.")}
            
            test <- sapply(flist,function(x) try(x(c(0,1,2,3,4,5)),silent=TRUE))
            if(!all(test %in% c("TRUE","FALSE"))) stop("All function in flist must take a vector of abundance values and return a logical.")
            
            criteria <- sapply(flist,function(x) apply(object@data,2,x))
            criteria <- as.logical(apply(criteria,1,prod))
            
            new.data <- object@data[,which(criteria),drop=F]
            new.taxa <- taxa(object)[which(criteria),,drop=F]
            
            if(join.trim){
              new.data <- cbind(new.data,"trim"=rowSums(object@data[,which(!criteria)]))
              new.taxa <- rbind(new.taxa,"trim"=rep("motley",nrank(object)))
            }
            
            return(mg(data=new.data,meta=meta(object),taxa=new.taxa))
          })
################################################################################
################################################################################
# END FILTER TAXA
################################################################################
################################################################################
