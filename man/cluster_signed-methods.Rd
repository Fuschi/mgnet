% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_signed.R
\docType{methods}
\name{cluster_signed}
\alias{cluster_signed}
\alias{cluster_signed,igraph-method}
\alias{cluster_signed,mgnet-method}
\alias{cluster_signed,list-method}
\title{Get Communities of a Signed Weighted Graph}
\usage{
cluster_signed(
  obj,
  OS = "Linux",
  Resistance = 0,
  Penalty_Coefficien = 1,
  add.names = TRUE
)

\S4method{cluster_signed}{igraph}(
  obj,
  OS = "Linux",
  Resistance = 0,
  Penalty_Coefficien = 1,
  add.names = TRUE
)

\S4method{cluster_signed}{mgnet}(obj, OS = "Linux", add.names = TRUE)

\S4method{cluster_signed}{list}(obj, OS = "Linux", add.names = TRUE)
}
\arguments{
\item{obj}{weighted undirected network belong to \code{\link{igraph}} class
or an mgnet.}

\item{OS}{(default Linux) string with the operating system running. Possible choices are
"Linux","Windows","Mac".}

\item{Resistance}{resistance of nodes to join communities, as a common
self-loop positive or negative real number, default no resistance with value set
to 0.}

\item{Penalty_Coefficien}{relative importance of null-case term non-negative
real number default set to 1.}

\item{add.names}{logical with default value set to TRUE. It indicates whether
you want to use the name of the vertices also in the resulting communities or
simply numerical indexes.}
}
\value{
\code{\link{communities}} igraph object able to manage to communities
graph info. The unique difference from igraph routine is the treatment with
the isolated nodes. In this case all isolated nodes are classified in the
community \code{'0'} and not as different communities of size one.
}
\description{
Adaptation in R of the algorithm for detention of weighted
communities with sign of a graph developed by Sergio Gomez.
https://deim.urv.cat/~sergio.gomez/radatools.php.
The algorithm tries to find dense subgraph, also called communities via
optimization of a signed definition of modularity score.
}
