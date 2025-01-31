% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2g_classes_gather.R
\name{gather_link}
\alias{gather_link}
\alias{gather_link,mgnet-method}
\alias{gather_link,mgnetList-method}
\title{Retrieve Edge List with Metadata from mgnet Objects}
\usage{
gather_link(object, .suffix)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.
For \code{mgnet}, the method extracts the network's edges and merges taxa metadata
(if requested).
For \code{mgnetList}, the method applies the transformation to each \code{mgnet} in
the list and combines the results, adding a \code{mgnet} column for the source.}

\item{.suffix}{A character vector of length 2 providing suffixes to append
to the taxa metadata columns for the two nodes in each link. The two suffixes
must be distinct (e.g. \code{c("_from", "_to")}). If missing or \code{NULL},
no metadata join is performed.}
}
\value{
A dataframe (tibble) containing the edges plus node-level metadata
(if requested). For \code{mgnetList}, the result is a single tibble with
an extra \code{mgnet} column identifying each sub-object.
}
\description{
Extracts edge information from the network slot of an \code{mgnet} or \code{mgnetList} object,
optionally merging metadata from both vertices (nodes) into each link.
Internally uses \code{igraph::as_data_frame()} to get the edge list, then
joins with taxa metadata if \code{.suffix} is provided.
}
\details{
\itemize{
\item \strong{Edge Attributes}: All available attributes on the edges (e.g., weight) are
preserved.
\item \strong{Node Metadata}: If \code{.suffix} is given, columns from the \code{from} and
\code{to} nodes' metadata are appended with the respective suffixes, so they
do not overwrite each other or edge attributes.
\item \strong{Filtering}: If the \code{mgnet} object has selected links (via \code{select_link}),
only those links are kept.
}
}
