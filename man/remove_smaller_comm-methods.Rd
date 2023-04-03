% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mgnet.R
\docType{methods}
\name{remove_smaller_comm}
\alias{remove_smaller_comm}
\alias{remove_smaller_comm,mgnet,numeric-method}
\alias{remove_smaller_comm,list,numeric-method}
\title{Remove smaller communities}
\usage{
remove_smaller_comm(obj, size, trim = TRUE)

\S4method{remove_smaller_comm}{mgnet,numeric}(obj, size, trim = TRUE)

\S4method{remove_smaller_comm}{list,numeric}(obj, size, trim = TRUE)
}
\arguments{
\item{obj}{mgnet.}

\item{size}{integer indicates the vertex number threshold}

\item{trim}{(Optional, default TRUE) logical. If true, the function removes
all nodes not belonging to a community with size equal to size. If false the
filtered vertices are set as isolated.}
}
\description{
Allows you to remove communities based on the number of vertices.
}