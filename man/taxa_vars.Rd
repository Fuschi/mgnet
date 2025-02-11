% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2d_classes_getters.R
\name{taxa_vars}
\alias{taxa_vars}
\alias{taxa_vars,mgnet-method}
\alias{taxa_vars,mgnetList-method}
\title{Get Taxa Metadata Variables}
\usage{
taxa_vars(object, .fmt = "list")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{.fmt}{Character; specifies the output format when checking an \code{mgnetList}.
Accepted values are \code{"list"} for a list of sample metadata variables
in each \code{mgnet} and \code{"unique"} for a single array with the unique elemnts.
Default is \code{"list"}.}
}
\value{
For an \code{mgnet} object, a character vector of metadata variable names.
For an \code{mgnetList} object, a named list of character vectors, with each list item representing
the metadata variable names in the corresponding \code{mgnet} objects.
}
\description{
Retrieves the names of metadata variables available in the \code{taxa} slot of an \code{mgnet} object,
or for each \code{mgnet} object within an \code{mgnetList}.
}
\examples{
data(mg, package = "mgnet")
taxa_vars(mg)  

data(mgl, package = "mgnet")
taxa_vars(mgl, .fmt = "list)  
taxa_vars(mgl, .fmt = "unique) 

}
