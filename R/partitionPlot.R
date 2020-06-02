#' Make matrix for layout to partition plotting area
#'
#' This function proposes a matrix for use with \code{\link[graphics]{layout}} to arrange given number of plots to be placed on a page/plotting area.
#' In certain instances the proposed layout may accomodate slightly more plots, eg \code{nFig=5} can not be arranged in 2 or 3 columns without an empty last spot.
#' Portrait (vertival) or lanscape (horizontal) layout proportions can be chosen.  The user can also impose a given number of columns.
#' 
#' @param nFig (integer) number of figures to be arrages on single plotting surface (ie window or plotting device)
#' @param returnMatr (logical) will return matrix ready for use by  \code{\link[graphics]{layout}}; returns vector with nRow and nCol if \code{=FALSE}
#' @param horiz (logical) will priviledge horizontal layout if \code{TRUE}
#' @param figNcol (integer) optional number of columns
#' @param byrow (logical) toggle if output is in order of rows or columns (equivament to \code{\link[base]{matrix}} 
#' @param callFrom (character) allows easier tracking of messages produced
#' @return matrix for use with \code{layout} or (if \code{returnMatr=FALSE} numeric vector with number of segements in x- an y-axis) 
#' @seealso \code{\link[graphics]{layout}}
#' @examples
#' partitionPlot(5); partitionPlot(14,horiz=TRUE)
#' @export  
partitionPlot <- function(nFig,returnMatr=TRUE,horiz=TRUE,figNcol=NULL,byrow=TRUE,callFrom=NULL) {
  ## function to create layout-matrix for multiple separate plots
  ## to do : test more non-square layouts for (near-)perfect fit to 'nFig'
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="partitionPlot")
  if(is.null(figNcol)) figNcol <- ceiling(sqrt(nFig))
  figNrow <- ceiling(nFig/figNcol)
  figDim <- c(ceiling(nFig/figNcol),figNcol)
  if(figDim[2] < figDim[1] & horiz | figDim[2] > figDim[1] & !horiz) figDim <- figDim[2:1]
  out <- if(returnMatr) matrix(1:(figDim[2]*figDim[1]),ncol=figDim[2],byrow=byrow) else figDim
  out }
    
