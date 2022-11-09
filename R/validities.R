################################################################################
setValidity("ngs_data", function(object){
  
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
setValidity("sample_metadata", function(object){
  
  if(any(dim(object)==0)) return("\n data.frame must have non-zero dimensions.")
  if(is.null(rownames(object))) return("\n matrix must have the rows names where the samples IDs were saved.")
  if(is.null(colnames(object))) return("\n matrix must have the cols names where the experimental variables were saved")
  if(any(duplicated(rownames(object)))) return("\n find in matrix at least a duplicated row name / sample ID.")
  if(any(duplicated(colnames(object)))) return("\n find in matrix at least a duplicated col name / experimental variable.")

  TRUE
})
################################################################################
setValidity("taxonomy_table", function(object){
  
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
setValidity("mg", function(object){
  
  if(!is.null(object@data) && !is.null(object@meta)){
    if(!all(rownames(object@data)==rownames(object@meta))) return("\n rows names / sample IDs must be identical in data and meta slots")
  }
  
  if(!is.null(object@data) && !is.null(object@taxa)){
    if(!all(colnames(object@data)==colnames(object@taxa))) return("\n cols names / taxa IDs must be identical in data and taxa slots")
  }
  
  TRUE
})