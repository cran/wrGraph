% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLinReg.R
\name{plotLinReg}
\alias{plotLinReg}
\title{Plot linear regression and confidence interval of regression}
\usage{
plotLinReg(
  dat,
  indepVarLst = NULL,
  dependVar = NULL,
  cusTxt = NULL,
  regrLty = 1,
  regrLwd = 1,
  regrCol = 1,
  confInt = 0.95,
  confCol = NULL,
  xLab = NULL,
  yLab = NULL,
  xLim = NULL,
  yLim = NULL,
  tit = NULL,
  nSignif = 3,
  col = 1,
  pch = 1,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{dat}{(numeric, data.frame or list) main data to plot/inspect. If numeric 'dat' will be used as dependent variable (y-data)
together with numeric 'indepVarLst' (independent variable); if list, then list-elments \code{indepVarLst} and \code{dependVar} will be used; if matrix, the the 1st and 2nd colum will be used}

\item{indepVarLst}{(character) if 'dat' is list, this designes the list element with the explanatory or independent variable (ie the variable used for explaining, typically x-data)}

\item{dependVar}{(character) if 'dat' is list, this designes the list element with dependent variable (ie the variable to be explained, typically y-data) to test}

\item{cusTxt}{(character) optional custom text to display in subtitle (instead of p-value to H0: slope.regression=0)}

\item{regrLty}{(integer) line type for regression}

\item{regrLwd}{(integer) line width for regression}

\item{regrCol}{(integer) color of regression-line}

\item{confInt}{(numeric, between 0 and 1) the probabiity alpha for the regression interval, if \code{NULL} no confidence intervall will be plotted/calculated}

\item{confCol}{(character) (background) color for confidence-interval}

\item{xLab}{(character) optional custom x-label}

\item{yLab}{(character) optional custom y-label}

\item{xLim}{(numeric) custom limit for x-axis (see also \code{\link[graphics]{par}})}

\item{yLim}{(numeric) custom limit for y-axis (see also \code{\link[graphics]{par}})}

\item{tit}{(character) optional title}

\item{nSignif}{(integer) number of significant digits for regression parameters in subtitle of plot}

\item{col}{(integer or character) custom color for points (choose \code{NULL} for not plotting the actual data)}

\item{pch}{(integer or character) type of symbol for points (see also \code{\link[graphics]{par}})}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This functions simply plots (to the current graphical devce); an invisible list containing $data, $linRegr, $confInterval (if calculated) may be returned, too
}
\description{
This function provides help to display a series of bivariate points given in 'dat' (multiple data formats possible), to model a linear regression and plot the results.
Furthermore, a confidence interval to the regression may be added to the plot, regression parameters get be displayed.
}
\examples{
set.seed(2020); dat1 <- rep(1:6,each=2) +runif(12,0,1)
plotLinReg(dat1, gl(6,2))
## extract elements out of list :
li2 <- list(aa=gl(5,2), bb=dat1[1:10])
plotLinReg(li2, indepVarLst="aa", dependVar="bb")
}
\seealso{
\code{\link[wrMisc]{exclExtrValues}} for decision of potential outliers; \code{\link[graphics]{hist}}, \code{\link{vioplotW}}
}
