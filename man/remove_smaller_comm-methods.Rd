% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mgnet.R
\docType{methods}
\name{remove_smaller_comm}
\alias{remove_smaller_comm}
\alias{remove_smaller_comm,mgnet,numeric,logical-method}
\alias{remove_smaller_comm,mgnet,logical}
\alias{remove_smaller_comm,list,numeric,logical-method}
\alias{remove_smaller_comm-methods,list,numeric,logical}
\title{Remove smaller communities}
\usage{
remove_smaller_comm(obj, size, trim)

\S4method{remove_smaller_comm}{mgnet,numeric,logical}(obj, size, trim)

\S4method{remove_smaller_comm}{list,numeric,logical}(obj, size, trim)
}
\arguments{
\item{obj}{\code{\link{mgnet-class}}}

\item{size}{integer indicates the vertex number threshold}

\item{trim}{logical. If true, the function removes all nodes not belonging
to a community with size equal to size. If false the filtered vertices are
set as isolated.}
}
\description{
allows you to remove communities based on the number of
vertices.
}
