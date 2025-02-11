% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a_classes_base_methods.R
\name{has_sample}
\alias{has_sample}
\alias{has_sample,mgnet-method}
\alias{has_sample,mgnetList-method}
\title{Check if Samples Are Present}
\usage{
has_sample(object, .fmt = "list")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{.fmt}{Character; specifies the output format when checking an \code{mgnetList}.
Accepted values are \code{"list"} for a list of logical values, \code{"any"} for
a single logical indicating if any \code{mgnet} contains samples, and \code{"all"}
for a single logical indicating if all \code{mgnet} objects contain samples. Default is \code{"list"}.}
}
\value{
For an \code{mgnet} object, a logical indicating whether samples are present.
For an \code{mgnetList} object, the output depends on the value of \code{.fmt}.
If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"},
returns a single logical value indicating if any \code{mgnet} contains samples.
If \code{.fmt = "all"}, returns a single logical value indicating if all \code{mgnet}
objects contain samples.
}
\description{
Checks whether there are samples present in an \code{mgnet} object or in each
\code{mgnet} object within an \code{mgnetList}. For \code{mgnetList} objects, the output format can
be a list, or a single logical value indicating if any or all objects contain samples.
}
\examples{
data(mg, package = "mgnet")
has_sample(mg)  

data(mgl, package = "mgnet")
has_sample(mgl, .fmt = "list")  
has_sample(mgl, .fmt = "any")  
has_sample(mgl, .fmt = "all")  

}
