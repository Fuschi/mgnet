################################################################################
setValidity("ngs_data", function(object){
  
  if(any(dim(object)==0)) return("\n matrix must have non-zero dimensions.")
  if(!is.numeric(object@.Data)) return("\n matrix must be numeric")
  if(!all(object@.Data>=0))return("\n all matrix elements must be greater or equal to zero")
  
  TRUE
})
################################################################################
setValidity("sample_metadata", function(object){
  
  if(any(dim(object)==0)) return("\n data.frame must have non-zero dimensions.")

  TRUE
})