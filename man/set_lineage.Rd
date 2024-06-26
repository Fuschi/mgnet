% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3c_classes_set.R
\name{set_lineage}
\alias{set_lineage}
\alias{set_lineage,mgnet-method}
\alias{set_lineage,mgnetList-method}
\title{Set \code{lineage} Slot in \code{mgnet} Objects}
\usage{
set_lineage(object, value)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{value}{The new lineage data to be set, a numeric matrix for \code{mgnet} objects
or a list of numeric matrices for \code{mgnetList} objects.}
}
\value{
The modified \code{mgnet} or \code{mgnetList} object with the updated lineage data.
}
\description{
This function sets the lineage data for an \code{mgnet} object or each \code{mgnet} object
within an \code{mgnetList}.
}
