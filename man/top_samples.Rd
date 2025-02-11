% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3e_classes_top.R
\name{top_samples}
\alias{top_samples}
\alias{top_samples,mgnet-method}
\alias{top_samples,mgnetList-method}
\title{Return Ordered the Sample IDs Based on Aggregate Metrics on Abundances Matrices}
\usage{
top_samples(
  object,
  field,
  .var,
  .aggregation_fun = sum,
  .order_fun = sum,
  decreasing = TRUE
)
}
\arguments{
\item{object}{An \code{mgnet} object.}

\item{field}{Character string specifying the data field to use for computing the metric.
Acceptable values are \code{"abun"}, \code{"rela"}, or \code{"norm"}.
The field must not be empty in the \code{mgnet} object.}

\item{.var}{A character string specifying the column name from the \code{taxa} slot to be used
for grouping and aggregation. If this parameter is omitted, the function will return
the original abundance data without any aggregation.}

\item{.aggregation_fun}{A function to specify how the abundance data should be aggregated.
This can be any function that takes a numeric vector and returns a single number (e.g., sum, mean).
The default aggregation function is sum.}

\item{.order_fun}{Function used to compute the aggregate metric across the selected field. The default
is \code{sum}, but other functions like \code{mean} or \code{median} can be used to tailor the analysis.}

\item{decreasing}{Logical indicating if sorting should be in decreasing order. When \code{TRUE} (default),
samples with the highest metric values are considered top.}
}
\value{
For \code{mgnet} objects, a character vector of sample IDs representing the ordered samples.
For \code{mgnetList} objects, returns a list where each element is a character vector of
the ordered samples for each respective \code{mgnet}.
}
\description{
Return the ordered samples from an \code{mgnet} object based on an aggregate metric calculated from
specified data fields such as abundance, relative abundance, or normalized abundance, with
the chosen metric, which can be customized.
}
\note{
Ensure that the specified field is not empty in the \code{mgnet} object before calling this function.
Failing to do so will result in an error.
}
