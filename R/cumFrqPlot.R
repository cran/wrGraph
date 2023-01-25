#' Cumulative (or sorted) frequency plot (takes columns of 'dat' as separate series)
#'
#' Display data as sorted or cumulative frequency plot. This type of plot represents an alternative to plotting data as histograms.
#' Histograms are very universal and which are very intuitive. However,fine-tuning the bandwith (ie width of the bars) may be very delicate, 
#' fine resultion details may often remain hidden. 
#' One of the advanges of directly displaying all data-points is that subtile differences may be revealed easier, compared to calssical histograms. 
#' Furthermore, the plot prensented her offeres more options to display multiple series of data simultaneaously. 
#' Thus, this type of plot may be useful to compare eg results of data normalization.
#' Of course, with very large data-sets (eg > 3000 values) this gain of 'details' will be less important (compared to histograms) and will penalize speed.
#' In such cases the argument \code{thisResol} will get useful as it allows to reduce the resultion and introduce binning.
#' Alternatively for very large data-sets one may looking into density-plots or vioplots (eg \code{\link{vioplotW}}).
#' The argument \code{CVlimit} allows optionally excluding extreme values. 
#' If numeric (& > 2 columns), its value will be used \code{\link[wrMisc]{exclExtrValues}} to identify series with column-median > 'CVlimit'.
#' Of course, exclusion of extreme values should be done with great care, important features of the data may get lost.
#'
#' @param dat (matrix or data.frame) data to plot/inspect
#' @param cumSum (logical) for either plotting cumulates Sums (then \code{thisResol} for number of breaks) or (if \code{=FALSE}) simply sorted values -> max resolution
#' @param exclCol (integer) columlns to exclude
#' @param colNames (character) for alternative column/series names in display, as long as \code{displColNa=TRUE}
#' @param displColNa (logical) display column-names
#' @param tit (character) custom title
#' @param xLim (numeric) custom limit for x-axis (see also \code{\link[graphics]{par}})
#' @param yLim (numeric) custom limit for y-axis (see also \code{\link[graphics]{par}})
#' @param xLab (character) custom x-axis label
#' @param yLab (character) custom y-axis label
#' @param col (integer or character) custom colors
#' @param CVlimit (numeric) for the tag 'outlier column' (uses \code{\link[wrMisc]{exclExtrValues}})  identify & mark column with median row-CV > \code{CVlimit}
#' @param thisResol (integer) resolution \code{res} for binning large data-sets
#' @param supTxtAdj (numeric) parameter \code{adj} for supplemetal text
#' @param supTxtYOffs (numeric) supplemental offset for text on y axis
#' @param useLog (character) default="", otherwise for setting axis in log-scale "x", "y" or "xy"
#' @param silent (logical) suppress messages
#' @param debug (logical) additonal messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function plots to the current garphical device
#' @seealso \code{\link[graphics]{layout}}, \code{\link[wrMisc]{exclExtrValues}} for decision of potential outliers; \code{\link[graphics]{hist}}, \code{\link{vioplotW}}
#' @examples
#' set.seed(2017); dat0 <- matrix(rnorm(500), ncol=5, dimnames=list(NULL,1:5))
#' cumFrqPlot(dat0, tit="Sorted values")
#' cumFrqPlot(dat0, cumSum=TRUE, tit="Sum of sorted values")
#' @export  
cumFrqPlot <- function(dat, cumSum=FALSE, exclCol=NULL, colNames=NULL, displColNa=TRUE, tit=NULL, xLim=NULL, yLim=NULL,
  xLab=NULL, yLab=NULL, col=NULL, CVlimit=NULL, thisResol=NULL, supTxtAdj=0, supTxtYOffs=0, useLog="", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## cumulative frequency plot (takes columns of 'dat' as separate series)      
  argNa <- deparse(substitute(dat))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="cumFrqPlot")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
if(length(dim(dat)) >2) {
    if(!silent) message(fxNa," ('dat' has >2 dims) Using ONLY 1st and 2nd dimension of 'dat' for plotting !")
    dat <- dat[,,1] }
  if(!is.null(colNames)) if(length(colNames) != ncol(as.matrix(dat))) {
    warning(" number of items in 'colNames' doesn't match data provided"); colNames <- NULL }
  if(length(exclCol) >0) dat <- dat[,-1*exclCol]
  if(is.null(colNames)) colNames <- colnames(dat)
  dat <- wrMisc::.keepFiniteCol(as.matrix(dat), silent, callFrom=fxNa)
  if(!is.matrix(dat) & !is.data.frame(dat)) dat <- as.matrix(as.numeric(dat))
  ## optional removing of extreme values (outlyers)
  if(length(CVlimit) >0 & is.numeric(CVlimit)) {
   dat <- apply(dat, 2, wrMisc::exclExtrValues, result="val", CVlim=CVlimit, maxExcl=1, goodValues=FALSE, silent=silent, callFrom=fxNa)
  }       
  ## sort
  dat <- apply(dat, 2, sort, decreasing=FALSE,na.last=TRUE)   # put NAs to end
  repNA <- !is.finite(dat)
  nNA <- colSums(repNA)
  if(any(nNA ==nrow(dat)))  {
    if(!silent) message(fxNa," problem with columns ",wrMisc::pasteC(colnames(dat)[which(nNA ==nrow(dat))],quote="'")," , removing")
    dat <- dat[,which(nNA <nrow(dat))]}    
  ## median for display
  datMed <- apply(dat,2,stats::median,na.rm=TRUE)
  ## take slices if thisResol
  ##  taking slices before cumSum should go much faster on very big data, but some loss in precision since need to mulitply results for cumSum
  if(length(thisResol) >0 & is.numeric(thisResol)) {
    if(thisResol[1] > nrow(dat)/1.5) {                  # can't go better than full resolution ! Or, if < 50% stay full res
      if(!silent) message(fxNa," reducing resolution to '",thisResol,"' brings not sufficient/major again, remain at full resolution")
      thisResol <- NULL }}   
  if(length(thisResol) >0) if(is.numeric(thisResol)) {
    xx <- round(nrow(dat)/thisResol)
    xy <- (1:thisResol)*xx
    chXY <- xy > nrow(dat)
    if(any(chXY)) xy <- xy[which(!chXY)]
    dat <- dat[xy,]
    if(cumSum) dat <- dat*xx
  }
  ## cumsum (if needed)
  if(cumSum) {
    dat <- apply(dat,2,cumsum)  ## note that NAs
  }  
  ##prepare figure
  if(length(col) <1) col <- 1:ncol(dat)
  if(is.null(xLab)) xLab <- paste("value of ", if(cumSum) "cumulated" else "sorted"," per column")
  if(is.null(yLab)) yLab <- "fraction of data"
  if(is.null(xLim)) xLim <- range(dat,na.rm=TRUE,finite=TRUE)
  if(is.null(yLim)) yLim <- c(0,1)  #c(1,nrow(dat))/nrow(dat)  
  if(length(xLim) !=2) {if(!is.null(xLim)) warning("invalid entry for 'xLim', ignoring"); xLim <- NULL }
  if(length(yLim) !=2) {if(!is.null(yLim)) warning("invalid entry for 'yLim', ignoring"); yLim <- NULL }  
  ## plot empty
  msg1 <- paste0(fxNa,": Unknow argument content ('",useLog,"') for 'useLog'; resetting to default no log")
  if(length(useLog) !=1) { if(!silent) message(msg1); useLog <- ""}
  if(!(useLog %in% c("","x","y","xy"))) { warning(msg1); useLog <- ""}
  graphics::plot(c(1,nrow(dat)), c(range(dat, na.rm=TRUE)),
    xlim=xLim, ylim=yLim, las=1, main=tit, xlab=xLab,ylab=yLab,type="n",log=useLog)                             # plot empty frame
  ## add lines and legend-text 
  for(j in 1:ncol(dat)) {
    graphics::lines(dat[,j],(1:nrow(dat))/nrow(dat), col=col[j])
    graphics::mtext(paste("  med ",if(displColNa) colNames[j],"=",signif(datMed[j],3)),
      line=-2*j/3+supTxtYOffs,cex=0.7,col=col[j],adj=supTxtAdj)
  }
}
  
