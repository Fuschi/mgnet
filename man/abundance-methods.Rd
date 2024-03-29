% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mgnet.R
\docType{methods}
\name{abundance}
\alias{abundance}
\alias{abundance,mgnet,missing-method}
\alias{abundance,mgnet,character-method}
\alias{abundance,list,missing-method}
\alias{abundance,list,character-method}
\title{Get abundances at choosen rank.}
\usage{
abundance(object,rank)

\S4method{abundance}{mgnet,missing}(object)

\S4method{abundance}{mgnet,character}(object, rank)

\S4method{abundance}{list,missing}(object)

\S4method{abundance}{list,character}(object, rank)
}
\arguments{
\item{object}{mgnet or list.}

\item{rank}{(Optional) character with the taxonomic rank choosen.}
}
\description{
Retrieves the abundances of data at choosen taxonomy rank. Abundance at
higher rank are returned as sums of abundance of elements with the same
classification.
}
