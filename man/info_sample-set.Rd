% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3b_classes_setters.R
\name{info_sample<-}
\alias{info_sample<-}
\alias{info_sample<-,mgnet-method}
\alias{info_sample<-,mgnetList-method}
\title{Set Sample Metadata}
\usage{
info_sample(object) <- value
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{value}{The new sample metadata to be set.}
}
\value{
The \code{mgnet} or \code{mgnetList} object with updated sample metadata.
}
\description{
This setter function updates the sample metadata for \code{mgnet} objects and
each \code{mgnet} object within a \code{mgnetList}. The metadata must be provided as a dataframe
for \code{mgnet} objects. For \code{mgnetList} objects, the metadata should be a named list of dataframes
corresponding to each \code{mgnet} object within the list.
}
