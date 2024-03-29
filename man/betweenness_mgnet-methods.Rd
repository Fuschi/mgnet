% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-mgnet.R
\docType{methods}
\name{betweenness_mgnet}
\alias{betweenness_mgnet}
\alias{betweenness_mgnet,mgnet-method}
\alias{betweenness_mgnet,list-method}
\title{Inta Communities Betweenness Centrality}
\usage{
betweenness_mgnet(obj, normalize = "none")

\S4method{betweenness_mgnet}{mgnet}(obj, normalize = "none")

\S4method{betweenness_mgnet}{list}(obj, normalize = "none")
}
\arguments{
\item{obj}{mgnet class.}

\item{normalize}{allows to normalize the betweenness centralities scores.
"nVertices" if it will be divided respect the vertex number, "L1" if the sum
constrain is equal to 1 and "none" to leave invariate.}
}
\description{
Calculates the betweenness centralities of each vertices inside
its communities. For the scope the function take into account only positive
links.
}
