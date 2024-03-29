% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkForLegLoc.R
\name{checkForLegLoc}
\alias{checkForLegLoc}
\title{Find best place on plot for placing legend}
\usage{
checkForLegLoc(
  matr,
  sampleGrp = NULL,
  showLegend = TRUE,
  suplSpace = 4,
  testCorner = 1:4,
  silent = TRUE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{matr}{(matrix, list or data.frame) main data of plot}

\item{sampleGrp}{(character or factor) with this option the text to be displayed in the legend may be taken into consideration for its length}

\item{showLegend}{(logical or character) decide if \code{matr} should be checked for best location; if \code{showLegend} contains any of the standard legend-location designations (eg 'topleft') it will be used in the output}

\item{suplSpace}{(numeric) allows to consider extra room taken in legend by symbol and surrounding space, interpreted as n additional characters}

\item{testCorner}{(integer) which corners should be considered (1=left-top, 2=right-top, right-bottom, left-bottom)}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additonal messages for debugging}

\item{callFrom}{(character) allows easier tracking of messages produced}
}
\value{
list with $showL indicating if legend is desired and $loc for the proposition of the best location, $nConflicts gives the counts of conflicts
}
\description{
This function tries to find the best location for placing a legend of a bivariate plot, ie scatter-plot.
All 4 corners of the data to plot are inspected for the least occupation by data plotted while displaying the content of \code{sampleGrp}.
Alternatively, by setting the argument \code{showLegend} the user-defined legend will be returned
}
\examples{
dat1 <- matrix(c(1:5,1,1:5,5), ncol=2)
grp <- c("abc","efghijk")
(legLoc <- checkForLegLoc(dat1, grp)) 
plot(dat1, cex=3)
legend(legLoc$loc, legend=grp, text.col=2:3, pch=1, cex=0.8)
}
\seealso{
\code{\link[graphics]{legend}}
}
