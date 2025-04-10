% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vioplotW.R
\name{vioplotW}
\alias{vioplotW}
\title{Violin-Plots Version W}
\usage{
vioplotW(
  x,
  ...,
  finiteOnly = TRUE,
  removeEmpty = FALSE,
  halfViolin = FALSE,
  boxCol = "def",
  hh = NULL,
  xlim = NULL,
  ylim = NULL,
  nameSer = NULL,
  cexNameSer = NULL,
  horizontal = FALSE,
  col = "rainbow",
  border = "black",
  xlab = NULL,
  ylab = NULL,
  cexLab = NULL,
  cexAxis = NULL,
  lty = 1,
  pointCol = NULL,
  cexPt = NULL,
  tit = NULL,
  las = 1,
  lwd = 1,
  rectCol = "black",
  at = 0,
  add = FALSE,
  wex = NULL,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{x}{(matrix, list or data.frame) data to plot, or first series of data}

\item{...}{(numeric) additional sets of data to plot}

\item{finiteOnly}{(logical) eliminate non-finite elements to avoid potential errors (eg when encountering \code{NA})}

\item{removeEmpty}{(logical) omit empty series (or less than 4 finite numeric entries) of data from plot}

\item{halfViolin}{(logical or character) decide with \code{TRUE} or \code{FALSE} if full or only half of violins should be plotted, if "pairwise" always 2 data-sets will be plotted back-to-back}

\item{boxCol}{(character) decide if boxplot should be adde inside the violin, use "def" for default transparent grey}

\item{hh}{(numeric, length <4) smoothing parameter (standard deviation to kernel function, if omited anormal optimal smoothing parameter is used); equivalent to argument \code{h} in package  \href{https://CRAN.R-project.org/package=vioplot}{vioplot} ; see also \code{\link[sm]{sm.density}}}

\item{xlim}{(\code{NULL} or numeric, length=2) custom limit on x-axis, see also \code{\link[graphics]{par}}}

\item{ylim}{(\code{NULL} or numeric, length=2) custom limit on y-axis, see also \code{\link[graphics]{par}}}

\item{nameSer}{(character) custom label for data-sets or columns (length must match number of data-sets)}

\item{cexNameSer}{(numeric) size of individual data-series labels as cex-expansion factor (see also \code{\link[graphics]{par}})}

\item{horizontal}{(logical) orientation of plot}

\item{col}{(character or integer) custom colors or gradients like 'rainbow', 'grayscale', 'heat.colors', 'topo.colors', 'Spectral' or 'Paired',  or you may use colors made by the package \href{https://CRAN.R-project.org/package=colorRamps}{colorRamps}}

\item{border}{(character) custom color for figure border}

\item{xlab}{(character) custom x-axis label}

\item{ylab}{(character) custom y-axis label}

\item{cexLab}{(numeric) size of axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})}

\item{cexAxis}{(numeric) size of numeric y-axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})}

\item{lty}{(integer) line-type for linear regression line (see also \code{\link[graphics]{par}})}

\item{pointCol}{(character or numeric) display of median: color (defauly white)}

\item{cexPt}{(numeric) display of median : size of point as cex-expansion factor (see also \code{\link[graphics]{par}})}

\item{tit}{(character) custom title to figure}

\item{las}{(integer) orientation of axis labels (see also \code{\link[graphics]{par}})}

\item{lwd}{(integer) width of line(s) (see also \code{\link[graphics]{par}})}

\item{rectCol}{(character) color of rectangle}

\item{at}{(numeric) custom locoation of data-series names, ie the points at which tick-marks are to be drawn, will be passed to \code{\link[graphics]{axis}}, it's length ust match the number of data-sets}

\item{add}{(logical) add to existing plot if \code{TRUE}}

\item{wex}{(integer) relative expansion factor of the violin}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This function plots a figure (to the current graphical device)
}
\description{
This function allows generating \href{https://en.wikipedia.org/wiki/Violin_plot}{Violin plots}) using a variety of input formats and offers additional options for colors.
Main input may be multiple vectors, a matrix or list of multiple data-elements (entries may be of variable length),
individual colors for different sets of data or color-gradients can be specified, and the display of n per set of data was integtated
(based on an inspiration from the discussion 'Removing-NAs-from-dataframe-for-use-in-Vioplot' on the forum Nabble).
It is also possible to plot pairwise half-violins for easier pairwise-comparisons (using \code{halfViolin="pairwise"}).
Many arguments are kept similar to \href{https://CRAN.R-project.org/package=vioplot}{vioplot} (here, the package \code{vioplot} is not required/used).
}
\details{
The (relative) width of the density-profiles ('Violins') may be manually adjusted using the parameter \code{wex} which applies to all profiles drawn.

Note : Arguments have to be given with full names, lazy evaluation of arguments will not work properly with this function (since '...' is used to capture additional data-sets).
Note : \href{https://CRAN.R-project.org/package=vioplot}{vioplot} offers better options for plotting formulas
}
\examples{
set.seed(2013)
dat6 <- matrix(round(rnorm(300) +3, 1), ncol=6,
	 dimnames=list(paste0("li",1:50), letters[19:24]))
vioplotW(dat6)
## variable number of elements (each n is displayed)
dat6b <- apply(dat6, 2, function(x) x[which(x < 5)])
dat6b[[4]] <- dat6b[[4]][dat6b[[4]] < 4]
vioplotW(dat6b, col="Spectral")
vioplotW(dat6b, col="Spectral" ,halfViolin="pairwise", horizontal=TRUE)
vioplotW(dat6b, col="Spectral", halfViolin="pairwise", horizontal=FALSE)
}
\seealso{
the package \href{https://CRAN.R-project.org/package=vioplot}{vioplot}, \code{\link[sm]{sm}} is used for the density estimation
}
