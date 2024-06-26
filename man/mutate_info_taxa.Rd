% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4e_classes_mutate.R
\name{mutate_info_taxa}
\alias{mutate_info_taxa}
\alias{mutate_info_taxa,mgnet-method}
\alias{mutate_info_taxa,mgnetList-method}
\title{Mutate Columns in info_taxa of mgnet Objects}
\usage{
mutate_info_taxa(object, ..., .by = NULL)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object to be modified.}

\item{...}{Transformations to apply, passed to dplyr::mutate().}

\item{.by}{Optional; a character vector specifying the columns in the \code{info_sample}
data frame on which to group the data before applying the filter conditions.
This allows for complex, group-based filtering operations such as filtering
within subsets defined by one or more attributes.}
}
\value{
The \code{mgnet} or \code{mgnetList} object with its \code{info_taxa} slot updated.
}
\description{
\code{mutate_info_taxa} leverages dplyr's mutate functionality to enable column transformations or additions within the \code{info_taxa}
data frame in \code{mgnet} objects. This can be particularly useful for adding new taxonomic information or adjusting existing ones.
}
\details{
This function allows users to modify existing columns or add new columns to the \code{info_taxa} slot of \code{mgnet} or \code{mgnetList} objects.
It utilizes dplyr's mutate semantics for flexible data transformation.
}
\seealso{
\code{\link[dplyr]{mutate}} for details on transformation conditions.
}
