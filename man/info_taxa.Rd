% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3a_classes_getters.R
\name{info_taxa}
\alias{info_taxa}
\alias{info_taxa,mgnet-method}
\alias{info_taxa,mgnetList-method}
\title{Get Taxa Information}
\usage{
info_taxa(object, .fmt = "df")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{.fmt}{A character string specifying the output format of the result.
Possible choices are "df" for data.frame, and "tbl" for tibble.
When ".fmt" is set to "tbl", for \code{mgnet} objects, the row names of the abundance matrix
are moved into a new column named \code{sample_id}, aligning with the reserved keyword in
\code{mgnet}.}
}
\value{
The content of the \code{info_taxa} slot for \code{mgnet} or a list of such contents for \code{mgnetList}.
}
\description{
Retrieves the taxa information stored in the \code{info_taxa} slot of an \code{mgnet} object
or each \code{mgnet} object within an \code{mgnetList}, with the option to format the output as
data.frame or tibble.
}
