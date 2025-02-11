% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2b_classes_link.R
\name{get_selected_links}
\alias{get_selected_links}
\alias{get_selected_links,mgnet-method}
\alias{get_selected_links,mgnetList-method}
\title{Get Selected Links from an mgnet Object}
\usage{
get_selected_links(object)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}
}
\value{
A numeric vector (or named list of numeric vectors, for \code{mgnetList})
containing the indices of selected links.
}
\description{
Retrieves the numeric indices of links previously selected by \code{\link[=select_link]{select_link()}}.
}
