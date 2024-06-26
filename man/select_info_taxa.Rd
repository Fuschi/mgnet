% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4d_classes_select_info.R
\name{select_info_taxa}
\alias{select_info_taxa}
\alias{select_info_taxa,mgnet-method}
\alias{select_info_taxa,mgnetList-method}
\title{Select Columns from info_taxa in mgnet Objects}
\usage{
select_info_taxa(object, ...)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object from which columns will be selected.}

\item{...}{Conditions specifying which columns to select, passed to dplyr::select().
This can include a variety of selectors like column names, dplyr helper functions (starts_with, ends_with, contains, etc.), or indices.
For detailed usage, refer to \code{\link[dplyr]{select}}.}
}
\value{
An \code{mgnet} object with the \code{info_taxa} slot updated to include only the selected columns, or an \code{mgnetList} object
where each contained \code{mgnet} object has its \code{info_taxa} slot similarly updated.
}
\description{
\code{select_info_taxa} uses dplyr's select functionality to enable precise selection of columns from the \code{info_taxa}
data frame in \code{mgnet} objects. This can be particularly useful for simplifying metadata before further analysis
or for focusing on a subset of metadata attributes.
}
\details{
This function allows users to select specific columns from the \code{info_taxa} slot of \code{mgnet} or \code{mgnetList} objects.
It utilizes dplyr's select semantics to provide a flexible interface for column selection based on column names or conditions.
}
\seealso{
\code{\link[dplyr]{select}} for details on the selection conditions.
\code{\link{mgnet}} and \code{\link{mgnetList}} for details on the object structures.
}
