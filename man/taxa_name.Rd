% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a_classes_base_methods.R
\name{taxa_name}
\alias{taxa_name}
\alias{taxa_name,mgnet,missing-method}
\alias{taxa_name,mgnet,character-method}
\alias{taxa_name,mgnetList,missing-method}
\alias{taxa_name,mgnetList,character-method}
\title{Retrieve Taxa Names at Specified Rank}
\usage{
taxa_name(object, rank = "missing")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{rank}{The taxonomic rank for which names are to be retrieved. If not specified,
the function returns names at the finest available taxonomic rank. For \code{mgnet} objects,
this is derived from the column names of the \code{lineage} matrix. For \code{mgnetList} objects,
it iterates over each contained \code{mgnet} object to perform the extraction.}
}
\value{
For an \code{mgnet} object, a character vector of taxa names at the specified rank.
For an \code{mgnetList} object, a list where each element is a character vector
of taxa names at the specified rank from each \code{mgnet} object within the list.
}
\description{
Extracts the names of taxa at the specified taxonomic rank from an \code{mgnet} object,
or for each \code{mgnet} object within an \code{mgnetList}. Unlike \code{taxa_id} which returns
identifiers associated with the row/column names of the input matrix/data.frame
(often reflecting the finest classification, such as OTUs, where IDs and names may coincide),
\code{taxa_name} delves into the \code{lineage} matrix to fetch descriptive names for the entities at the
chosen rank. This provides more characteristic names of the entities.
}
