% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_reshape.R
\name{long_abundance_join}
\alias{long_abundance_join}
\title{Join and Convert Abundance Data to Long Format}
\usage{
long_abundance_join(object, needed_keys)
}
\arguments{
\item{object}{An \code{mgnet} object.}

\item{needed_keys}{A list containing the abundance keys to be used (e.g., \code{abun}, \code{rela}, \code{norm}).}
}
\value{
A tibble with \code{sample_id}, \code{taxa_id}, and the abundance variables in long format.
}
\description{
This internal function creates a grid of \code{sample_id} and \code{taxa_id} combinations
and joins all required abundance variables in long format. It converts abundance data
(e.g., \code{abun}, \code{rela}, \code{norm}) from wide to long format for easier manipulation.
}
\keyword{internal}
