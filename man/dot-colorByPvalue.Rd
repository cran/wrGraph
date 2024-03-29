% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VolcanoPlotW.R
\name{.colorByPvalue}
\alias{.colorByPvalue}
\title{Colors based on p-Values}
\usage{
.colorByPvalue(
  x,
  br = NULL,
  col = NULL,
  asIndex = FALSE,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{x}{(numeric) p-values (main input)}

\item{br}{(numeric) custom breaks (used with cut)}

\item{col}{custom colors (must be of length(br) -1)}

\item{asIndex}{(logical) custom breaks (used with cut)}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) supplemental messages for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This function retruns a color-gradient based on p-Values
}
\description{
This function helps defining color-gradient based on p-Values.
This fuction requires package RColorBrewer being installed
}
\examples{
.colorByPvalue((1:10)/10)
}
\seealso{
(for PCA) \code{\link{plotPCAw}})
}
