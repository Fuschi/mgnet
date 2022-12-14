% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mg.R
\docType{methods}
\name{selection_taxa}
\alias{selection_taxa}
\alias{selection_taxa,mg,vector-method}
\alias{selection_taxa,mg,vector}
\alias{selection_taxa,list,vector-method}
\alias{selection_taxa,list,vector}
\title{Selection a subset of taxa from mg object.}
\usage{
selection_taxa(object, idx)

\S4method{selection_taxa}{mg,vector}(object, idx)

\S4method{selection_taxa}{list,vector}(object, idx)
}
\arguments{
\item{object}{(Required) \code{\link{mg-class}}.}

\item{idx}{(Required) Vector of integer position indices or logical or names
of taxa (like to \code{[} extractor function).}
}
\description{
The function takes as input a vector of logical or position
indices to evaluate which taxa to keep.
}
