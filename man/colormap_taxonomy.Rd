% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility-color.R
\name{colormap_taxonomy}
\alias{colormap_taxonomy}
\title{Create color palette for taxonomy}
\usage{
colormap_taxonomy(
  taxaID,
  distinctColor = 20,
  extraColor = rgb(0.75, 0.75, 0.75, 0.25),
  alpha = 1,
  colorspace = "pretty"
)
}
\arguments{
\item{taxaID}{character vector with all taxonomic classification.}

\item{distinctColor}{max number of different final colors. After the values imposed
with distinctColor the function associated a light grey to every taxaID.}

\item{extraColor}{color of out of range colors, by default is light gray.}

\item{alpha}{(Optional) Level of trasparency of colors.}

\item{colorspace}{(Optional) Permits to set the color ranges. For more details
see \code{\link{qualpal}}.}
}
\description{
User wrapper for distinctColorPalette function of randomcoloR
package. Create named vector with taxaID as name and colors as values. If are
present duplicated in the vector the function choice only unique values.
}
