% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkForLegLoc.R
\name{.bestLegendLoc}
\alias{.bestLegendLoc}
\title{Search best corner of plot for placing for legend 
 
This function aims to find best corner for plotting a legend.}
\usage{
.bestLegendLoc(
  dat,
  txtLen,
  txtHi,
  displayPlotSearch = FALSE,
  silent = TRUE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{dat}{(matrix, list or data.frame) main data of plot}

\item{txtLen}{(numeric, length=1)}

\item{txtHi}{(numeric, length=1) text height from graphics::strheight() (including inter-line)}

\item{displayPlotSearch}{(logical) decide if lines to mark area where data is searched for legend should be drawn}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additonal messages for debugging}

\item{callFrom}{(character) allows easier tracking of messages produced}
}
\value{
numeric vector with counts of umber of points expected to enter legend-location for each corner (ie legend-localization)
}
\description{
Search best corner of plot for placing for legend 
 
This function aims to find best corner for plotting a legend.
}
\examples{
dat1 <- matrix(c(1:5,1,1:5,5), ncol=2)
(legLoc <- .bestLegendLoc(dat1, txtLen=0.4, txtHi=28))
}
\seealso{
\code{\link{checkForLegLoc}},  \code{\link[graphics]{legend}}
}
