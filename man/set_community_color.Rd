% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7c_classes_taxa_color.R
\name{set_community_color}
\alias{set_community_color}
\alias{set_community_color,mgnet-method}
\alias{set_community_color,mgnetList-method}
\title{Set Color Mapping to Community in an \code{mgnet} Object or \code{mgnetList}}
\usage{
set_community_color(
  object,
  sizes_threshold,
  smaller_color = "#D3D3D3",
  isolated_color = "#FFFFFF",
  alpha = 1,
  colorspace = "pretty",
  color_to,
  distinct_color = TRUE
)
}
\arguments{
\item{object}{An \code{mgnet} or \code{mgnetList} object.}

\item{sizes_threshold}{An integer specifying the minimum size a community must have to be assigned a unique color.
Communities with sizes below this threshold will be colored using \code{smaller_color}.}

\item{smaller_color}{A hexadecimal color code representing the color to be used for smaller communities, i.e.,
those communities whose sizes do not meet the \code{sizes_threshold}. Defaults to light gray ("#D3D3D3").}

\item{isolated_color}{A hexadecimal color code used to represent isolated nodes, typically assigned to community '0' in mgnet.
Defaults to white ("#FFFFFF"). See \link[mgnet]{as_mgnet_communities}.}

\item{alpha}{A numeric value between 0 and 1 that sets the transparency level of the colors, where 1 is fully opaque and
0 is fully transparent.}

\item{colorspace}{Defines the scheme used to generate the color palette. Possible options are 'pretty', 'pretty_dark',
'rainbow', 'pastels', or a custom list specifying the 'h' (hue), 's' (saturation), and 'l' (lightness) components.
Each option provides a different aesthetic and can be chosen based on the visual requirements of the analysis.
The function encapsulate qualpal to generate qualitative distinct color, see \link[qualpalr]{qualpal} for more details.}

\item{color_to}{The name of the new column within \code{info_taxa} where the assigned colors will be stored.}

\item{distinct_color}{Logical indicating whether distinct colors should be used for each \code{mgnet} object within an \code{mgnetList}
when \code{object} is an \code{mgnetList}. If TRUE, each \code{mgnet} object's communities are assigned unique colors independently;
if FALSE, a shared palette is used for all communities across the \code{mgnetList}, based on the maximum community number.}
}
\value{
An updated \code{mgnet} or \code{mgnetList} object with the \code{info_taxa} slot modified to include a new column named as specified
by \code{color_to}, which contains the color codes assigned based on community membership and size criteria.
}
\description{
This function assigns a visually distinct color palette to the communities detected in an \code{mgnet} object
or across multiple \code{mgnet} objects within an \code{mgnetList}. It provides options to handle communities based on their sizes,
assigning unique colors to larger communities while using a common color for smaller ones, thus facilitating clearer
visual distinctions in subsequent analyses or visualizations.
}
\details{
The function first assesses the community structure within the provided \code{mgnet} or \code{mgnetList} object. It applies a color
generation process that respects the specified \code{sizes_threshold}, assigning distinct colors to larger communities and a common
color to smaller ones. The \code{isolated_color} is used specifically for nodes that do not belong to any detected community or are
considered outliers. This function is particularly useful in network analyses where community detection results need to be
visualized to understand the structure and composition of the network.
}
\seealso{
\link[mgnet]{colormap_communities} \link[mgnet]{cluster_signed}
}
