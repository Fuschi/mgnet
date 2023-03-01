#' #' Spiec-Easi MB weights extraction
#' #' 
#' #' @description Extract the weight matrix of the Meinshausen and Bühlmann criterion
#' #' in the Spiec-Easi algorithm.
#' #' @importFrom SpiecEasi symBeta getOptBeta
#' #'
#' #' @param res.mb Spiec-Easi algorithm results with mb criteria
#' #'
#' #' @export
#' get_mb_weights <- function(res.mb){
#'   
#'   #CHECK ARGUMENTS
#'   if(class(res.mb)!="pulsar.refit") stop("class must be pulsar.refit")
#'   if(res.mb$est$method!="mb") stop("must be mb method")
#'   #END CHECK
#'   
#'   # get weights matrix W
#'   W <- as.matrix(symBeta(as.matrix(getOptBeta(res.mb)), mode='maxabs'))
#'   colnames(res.mb$est$data) -> colnames(W) -> rownames(W)
#'   W[W>0] <- 1
#'   W[W<0] <- -1
#'   
#'   return(W)
#' }
#' 
#' #' Spiec-Easi GLASSO weights extraction
#' #' 
#' #' @description Extract the weight matrix of the GLASSO criterion
#' #' in the Spiec-Easi algorithm.
#' #' 
#' #' @importFrom SpiecEasi getOptCov getRefit
#' #' @importFrom stats cov2cor
#' #'
#' #'@export
#' get_gl_weights <- function(res.gl, Refit=TRUE){
#'   
#'   #CHECK ARGUMENTS
#'   if(class(res.gl)!="pulsar.refit") stop("class must be pulsar.refit")
#'   if(res.gl$est$method!="glasso") stop("must be gl method")
#'   #END CHECK
#'   
#'   # get weights matrix W
#'   W <- cov2cor(as.matrix(getOptCov(res.gl)))
#'   if(Refit){W <- W * as.matrix(getRefit(res.gl))}
#'   
#'   colnames(res.gl$est$data) -> colnames(W) -> rownames(W)
#'   
#'   return(W)
#' }