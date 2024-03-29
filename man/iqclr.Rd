% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clr-transformations.R
\name{iqclr}
\alias{iqclr}
\title{Robust Centered Log-Ratio Transformation}
\usage{
iqclr(X, mar = 1)
}
\arguments{
\item{X}{(Required) numeric matrix or vector with all elements greater than or equal to 0.}

\item{mar}{(Optional) Integer giving the dimension where the function will be applied;
1 for rows and 2 for columns (default 1).}
}
\description{
It calculates the interquantile centered log-ratio transformation of X.
}
