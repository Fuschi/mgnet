% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a_classes_base_methods.R
\name{sample_id}
\alias{sample_id}
\alias{sample_id,mgnet-method}
\alias{sample_id,mgnetList-method}
\title{Get sample IDs from mgnet or mgnetList Objects}
\usage{
sample_id(object, .fmt = "list")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object from which to extract sample IDs.}
}
\value{
Retrieves the IDs of sample from an \code{mgnet} object or lists of sample IDs for each
\code{mgnet} object within an \code{mgnetList}.
}
\description{
Get sample IDs from mgnet or mgnetList Objects
}
\examples{
data(mg, package = "mgnet")
sample_id(mg)  

data(mgl, package = "mgnet")
sample_id(mgl)  

}
\seealso{
\link{mgnet} and \link{mgnetList} for details on the classes.
}
