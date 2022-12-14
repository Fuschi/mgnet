% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mg.R
\docType{methods}
\name{selection_sample}
\alias{selection_sample}
\alias{selection_sample,mg,vector-method}
\alias{selection_sample,mg,vector}
\alias{selection_sample,list,vector-method}
\alias{selection_sample,list,vector}
\title{Selection a subset of samples from mg object.}
\usage{
selection_sample(object, idx)

\S4method{selection_sample}{mg,vector}(object, idx)

\S4method{selection_sample}{list,vector}(object, idx)
}
\arguments{
\item{object}{(Required) \code{\link{mg-class}}.}

\item{idx}{(Required) Vector of integer position indices or logical (like to
\code{[} extractor function).}
}
\description{
The function takes as input a vector of logical or position
indices to evaluate which samples to keep.
}
