% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mg.R
\docType{methods}
\name{aggregate_taxa}
\alias{aggregate_taxa}
\alias{aggregate_taxa,mg,character-method}
\alias{aggregate_taxa,mg,character}
\alias{aggregate_taxa,list,character-method}
\alias{aggregate_taxa,list,character}
\title{Organize data in higher taxonomic level.}
\usage{
aggregate_taxa(object, rank)

\S4method{aggregate_taxa}{mg,character}(object, rank)

\S4method{aggregate_taxa}{list,character}(object, rank)
}
\arguments{
\item{object}{(Required) \code{\link{mg-class}}.}

\item{rank}{taxonomic level choosen.}
}
\description{
Reorganize an \code{\link{mg-class}} object in an higher taxonomy rank.
The function sums the taxa with the same classification and return a new mg
object with appropriate data and taxa slots.
}
