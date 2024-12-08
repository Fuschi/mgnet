% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a_classes_base_methods.R
\name{miss_slot}
\alias{miss_slot}
\alias{miss_slot,mgnet-method}
\alias{miss_slot,mgnetList-method}
\title{Check if a Slot is Missing in \code{mgnet} or \code{mgnetList} Objects}
\usage{
miss_slot(object, slot_name, .fmt = "list")
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{slot_name}{A character string specifying the slot to check for absence
(e.g., \code{"meta"}, \code{"taxa"}, \code{"netw"}, \code{"comm"}, \code{"abun"}, \code{"rela"}, \code{"norm"}).}

\item{.fmt}{Character; specifies the output format when checking an \code{mgnetList}.
Accepted values are \code{"list"} for a list of logical values, \code{"any"} for
a single logical indicating if any \code{mgnet} is missing the slot, and \code{"all"}
for a single logical indicating if all \code{mgnet} objects are missing the slot.
Default is \code{"list"}.}
}
\value{
For an \code{mgnet} object, a logical indicating whether the slot is missing.
For an \code{mgnetList} object, the output depends on the value of \code{.fmt}.
If \code{.fmt = "list"}, returns a named logical vector. If \code{.fmt = "any"},
returns a single logical value indicating if any \code{mgnet} is missing the slot.
If \code{.fmt = "all"}, returns a single logical value indicating if all \code{mgnet}
objects are missing the slot.
}
\description{
This function checks whether a specific slot (e.g., \code{meta}, \code{taxa}, \code{netw}, etc.)
is missing (i.e., empty or absent) in an \code{mgnet} object or in each \code{mgnet} object within an \code{mgnetList}.
For \code{mgnetList} objects, the output format can be a list, or a single logical
value indicating if any or all objects are missing the specified slot.
}
\examples{
data(mg, package = "mgnet")
miss_slot(mg, "meta")  

data(mgl, package = "mgnet")
miss_slot(mgl, "taxa", .fmt = "list")  
miss_slot(mgl, "netw", .fmt = "any")  
miss_slot(mgl, "comm", .fmt = "all")  

}
