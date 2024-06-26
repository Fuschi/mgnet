% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/executables.R
\name{signed-community-detection-executables}
\alias{signed-community-detection-executables}
\title{Executable Files for Signed Community Detection}
\description{
This package includes executable files that implement the algorithm for signed community detection developed by Sergio Gomez. These executables are used to analyze networks and identify communities within them, considering both the weight and sign of the connections between nodes.
}
\details{
The executables are utilized by specific functions within this package to perform community detection on networks constructed from metagenomic data. The results of this analysis can help in understanding the structure and dynamics of microbial communities.
}
\section{Executables Included}{

The package contains the following executables for different operating systems:
\itemize{
\item \code{Communities_Detection_Linux.exe} for Linux.
\item \code{Communities_Detection_Windows.exe} for Windows.
\item \code{Communities_Detection_Mac.exe} for macOS.
}
}

\section{Usage}{

Users do not need to interact directly with these executables. Instead, they are called automatically by the package's functions when performing community detection analyses. Ensure that the package is properly installed, and these executables will be utilized as needed.
}

\section{Source and License}{

These executables are based on the algorithm developed by Sergio Gomez. For more information about the algorithm and its implementation, visit \href{https://deim.urv.cat/~sergio.gomez/radatools.php}{Sergio Gomez's website}.

The use of these executables within this package complies with their licensing terms. Users are encouraged to refer to the original source for detailed licensing information and to ensure compliance with any usage restrictions.
}

