% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4c_classes_filter.R
\name{filter_taxa}
\alias{filter_taxa}
\alias{filter_taxa,mgnet-method}
\alias{filter_taxa,mgnetList-method}
\title{Filter \code{mgnet} or \code{mgnetList} Objects by Taxa Metadata and Abundance Fields}
\usage{
filter_taxa(object, ..., .by)

\S4method{filter_taxa}{mgnet}(object, ..., .by)

\S4method{filter_taxa}{mgnetList}(object, ..., .by)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object. For \code{mgnetList}, each contained
\code{mgnet} can be filtered according to the combined conditions, with a column
\code{mgnet} distinguishing them in any joined data frames.}

\item{...}{Unquoted expressions for filtering. These expressions can reference
columns in the taxa metadata (\code{meta}) or any of the abundance-related fields
(\code{abun}, \code{rela}, \code{norm}). Each expression is applied sequentially, and only
taxa that satisfy \strong{all} expressions are retained.}

\item{.by}{A character vector of grouping columns. The default differs by object
type:
\itemize{
\item For \code{mgnet}, the default is \code{.by = "taxa_id"}.
\item For \code{mgnetList}, the default is \code{.by = c("mgnet", "taxa_id")}.
}
If \code{.by} is \code{NULL}, no grouping is applied, and conditions are evaluated
across all taxa globally. If \code{.by} is missing, any existing taxa grouping
(from \code{group_taxa()}) is used instead.}
}
\value{
An \code{mgnet} or \code{mgnetList} object containing only the taxa (and their
corresponding abundance data, if any) that satisfy all provided filter conditions.
The rest are removed. If no taxa match, an object with zero taxa is returned
(but the same structure otherwise).
}
\description{
Filters the taxa in an \code{mgnet} or \code{mgnetList} object based on user-specified
conditions provided through \code{...}. Each condition can reference columns in the
taxa metadata (\code{taxa}) or abundance-related fields (\code{abun}, \code{rela}, \code{norm}).
When multiple conditions are provided, they are combined using a logical AND,
meaning a taxon must satisfy \strong{all} conditions to be retained.
}
\details{
This function integrates closely with \strong{tidyverse} conventions. Each expression
in \code{...} is captured and evaluated using tidy evaluation, allowing you to write
conditions similarly to \code{dplyr::filter()}. When an expression references
abundance variables (\code{abun}, \code{rela}, or \code{norm}), the function automatically
joins the taxa metadata with the relevant abundance data (\code{long_abun}) before
filtering.

\strong{Grouping Logic}:
\enumerate{
\item \strong{If you provide an explicit \code{.by}} (a non-empty character vector), those
columns define the grouping.
\item \strong{If \code{.by} is \code{NULL}}, no grouping is applied.
\item \strong{If \code{.by} is missing}:
\itemize{
\item \strong{If the object has existing taxa grouping} (set via \code{group_taxa()}), that
grouping is used.
\item \strong{Else}, if any filter expression references abundance variables,
grouping defaults to:
\itemize{
\item \code{c("mgnet", "taxa_id")} for \code{mgnetList} objects.
\item \code{"taxa_id"} for \code{mgnet} objects.
}
\item \strong{Otherwise}, no grouping is applied.
}
}

\strong{Multiple Conditions}:
Each condition in \code{...} is applied sequentially, narrowing down the set of
remaining taxa. After processing all conditions, only those taxa that
satisfied \strong{all} conditions remain.

\strong{Behavior with \code{mgnetList}}:
\itemize{
\item The \code{mgnet} column in the joined data frames identifies which \code{mgnet} each taxon
belongs to.
\item Filtering expressions can apply to metadata or abundance columns, referencing
each \code{mgnet}'s data.
\item Taxa are filtered \strong{independently} for each \code{mgnet} based on the conditions.
}

\strong{Return Value}:
\itemize{
\item For an \code{mgnet} object, an updated \code{mgnet} containing only the taxa that
passed all filters.
\item For a \code{mgnetList}, an updated list of \code{mgnet} objects. Each \code{mgnet} is
subset to the remaining taxa.
}

This approach provides a flexible, \code{dplyr}-like interface for complex filtering
of both taxa metadata and abundance data, respecting the specified grouping context.
}
