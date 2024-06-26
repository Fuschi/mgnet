% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2c_classes_extractors.R
\name{[,mgnet-method}
\alias{[,mgnet-method}
\title{Subset mgnet Object}
\usage{
\S4method{[}{mgnet}(x, i, j)
}
\arguments{
\item{x}{An \code{mgnet} object to be subsetted.}

\item{i}{Indices or logical vector indicating samples to be included in the subset. If missing, all samples are included.}

\item{j}{Indices or logical vector indicating taxa to be included in the subset. If missing, all taxa are included.}
}
\value{
A new \code{mgnet} object containing only the specified subset of samples and taxa. The structure and type of data within the \code{mgnet} object are preserved.
}
\description{
Subsets an \code{mgnet} object based on provided sample and taxa indices.
This method allows for the selective extraction of abundance, log abundance data,
sample and taxa metadata, and optionally, network and community data based on
specified indices.
}
\details{
The method enables the extraction of specific parts of the \code{mgnet}
object based on sample and taxa indices. When subsetting by samples, the abundance,
log abundance, and sample metadata are filtered accordingly.
Subsetting by taxa filters the taxa-related data and, if present, updates the
network and community slots to reflect the selected taxa.
It's important to note that subsetting by samples remove the network and community slots,
and a warning will be sended in such cases.
}
\seealso{
\code{\link[base]{Extract}} for base subsetting functions.
}
