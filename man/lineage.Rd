% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3a_classes_getters.R
\name{lineage}
\alias{lineage}
\alias{lineage,mgnet-method}
\alias{lineage,mgnetList-method}
\title{Get Lineage Information}
\usage{
lineage(object, .fmt = "mat")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{.fmt}{A character string specifying the output format of the result.
Possible choices are "mat" for matrix, df" for data.frame, and "tbl" for tibble.
When ".fmt" is set to "tbl", for \code{mgnet} objects, the row names of the abundance matrix
are moved into a new column named \code{sample_id}, aligning with the reserved keyword in
\code{mgnet}. Default is "mat".}
}
\value{
The content of the \code{lineage} slot for \code{mgnet} or a list of such contents for \code{mgnetList}.
}
\description{
Retrieves the taxonomic lineage information stored in the \code{lineage} slot of an \code{mgnet} object
or each \code{mgnet} object within an \code{mgnetList}, with the option to format the output as
matrix, data.frame or tibble.
}
