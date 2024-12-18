% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7a_classes_color.R
\name{set_community_color}
\alias{set_community_color}
\alias{set_community_color,mgnet-method}
\alias{set_community_color,mgnetList-method}
\title{Set Community Colors in \code{mgnet} Objects}
\usage{
set_community_color(
  object,
  size,
  color_to = "color_comm",
  smaller_color = "#D3D3D3",
  isolated_color = "#FFFFFF",
  colorspace = "pretty",
  alpha = 1
)
}
\arguments{
\item{object}{An \code{mgnet} object that includes community data.}

\item{size}{An integer that specifies the minimum size a community must have to receive a unique color.
Communities smaller than this size will be colored with \code{smaller_color} (default 1).}

\item{color_to}{The name of the column in the taxa data where the color values will be stored (default 'color_comm').}

\item{smaller_color}{A hexadecimal color code used for communities that do not exceed the \code{size} threshold,
defaulting to a light gray ("#D3D3D3").}

\item{isolated_color}{A hexadecimal color code used for isolated nodes, often those not belonging to any community,
defaulting to white ("#FFFFFF").}

\item{colorspace}{The qualitative color palette from which to generate the community colors. Options include
'pretty', 'pretty_dark', 'rainbow', 'pastels', and others as supported by the \code{qualpalr} package. This
parameter allows customization of the color scheme to fit different visualization needs.}

\item{alpha}{A numeric value between 0 and 1 indicating the opacity of the colors, where 1 is fully opaque and
0 is completely transparent. This is useful for creating layered visual effects in plots.}
}
\value{
The modified \code{mgnet} object with the added or updated community color data.
}
\description{
\code{set_community_color} assigns colors to taxa based on community memberships that are specified in the \code{comm_id}
of the \code{mgnet} object. The function employs a color palette that is sensitive to community size, differentiating
between larger and smaller communities.
}
\details{
This function applies a color scheme to taxa based on community membership, which is determined by
the sizes of communities identified in a metagenomic network. Each community can be assigned a unique
color, while smaller communities and isolated nodes are assigned default colors to maintain visual
clarity.

The function modifies the provided \code{mgnet} object by adding or updating a column with community colors.
It relies on the \code{colormap_community} function to generate a suitable color palette based on community
sizes and the specified colorspace. The visualization can then be tailored to highlight community structures
within metagenomic networks effectively.
}
\seealso{
\link[colormap_community]{Generate Color Palette for Community Visualization}
}
