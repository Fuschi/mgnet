% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4c_classes_filter.R
\name{filter_taxa}
\alias{filter_taxa}
\alias{filter_taxa,mgnet-method}
\alias{filter_taxa,mgnetList-method}
\title{Filter \code{mgnet} or \code{mgnetList} Objects Based on Taxa and Abundance Data}
\usage{
filter_taxa(object, ..., .by = NULL)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.
This function operates on the entire structure of these objects, allowing for the filtering
of taxa based on both metadata in the taxa field and various forms of abundance data.}

\item{...}{Dynamic expressions or functions to apply for filtering.
Users can specify conditions involving both metadata and abundance data, facilitating
complex and targeted filtering strategies that are essential in data-driven studies.}

\item{.by}{Optional; a character vector specifying the columns to group data by before transformations.
Defaults to 'taxa_id' for mgnet objects and c('mgnet', 'taxa_id') for mgnetList objects.
Grouping ensures that transformations are contextually applied within each subgroup defined
by .by. If you do not wish to group the data, set .by to NULL.}
}
\value{
Returns the \code{mgnet} or \code{mgnetList} object with updated content reflecting the applied filters.
This includes any transformations to taxa metadata and adjustments based on abundance data.
The structure outside the targeted data fields remains unchanged, preserving the integrity of the objects.
}
\description{
This function filters entire \code{mgnet} or \code{mgnetList} objects, including their taxa and abundance data,
based on user-specified conditions. It leverages the full suite of \code{tidyverse} tools, particularly
\code{dplyr}, to enable powerful and flexible data transformations.
}
\details{
Leveraging the \code{tidyverse} functionalities, this function allows for the integration of
sophisticated data manipulation techniques. It supports conditional operations, group-based
filtering, making it highly effective for analyses where detailed
taxa selection and data refinement are required.

\if{html}{\out{<div class="sourceCode">}}\preformatted{     ### Keywords in `mgnet` and `mgnetList`:
     - **abun, rela, norm**: Slots within `mgnet` objects that store abundance data, which can be
       directly manipulated or used in conjunction with metadata to perform advanced analyses.
     - **taxa_id**: An essential identifier used to uniquely reference individual taxa within an `mgnet` object. 
     - **mgnet**: Used exclusively within `mgnetList` objects to differentiate between multiple `mgnet` objects 
       contained in the list supporting sophisticated multi-dataset management.
}\if{html}{\out{</div>}}
}
