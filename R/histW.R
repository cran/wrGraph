#' Histogram (version by WR)
#'
#' This function proposes a few special tweaks to the general \code{\link[graphics]{hist}} function :
#' In a number of settings data are treated and plotted as log-data. This function allows feeding directly log2-data and displaying the x-axis   
#' (re-translated) in linear scale (see argument \code{isLog}). 
#' The default settings allow making (very) small histograms ('low resolution'), which may be used as a rough overview of bandwidth and distribution of values in \code{dat}. 
#' Similar to \code{\link[graphics]{hist}}, by changing the parameters \code{nBars} and/or \code{breaks} very 'high resolution' histograms can be produced.
#' By default it displays n per set of data (on the top of the figure). 
#' Note that the argument for (costom) title \code{main} is now called \code{tit}.
#' 
#' @param dat (matrix, list or data.frame) data to plot
#' @param fileName (character) name of file for saving graphics
#' @param output (character, length=1) options for output on 'screen' or saving image in various formats (set to 'jpg','png' or 'tif')
#' @param nBars (integer) number of bars in histogram (default for 'low resolution' plot to give rough overview)
#' @param breaks (integer) for (partial) compatibility with hist() : use only for number of breaks (or 'FD'), gets priority over 'nBars'
#' @param tit (character) custom title
#' @param subTi (character) may be \code{FALSE} for NOT displaying, or any text, otherwise range
#' @param xLab (character) custom x-axes label
#' @param yLab (character) custom y-axes label 
#' @param las (integer) optional fixed text orientation of x-axis numbers : use 1 for horizontal and 2 for perpendicular, see also \code{\link[graphics]{par}} 
#' @param xcex (numeric) cex-type expansion factor for x-axis numbers, see also \code{\link[graphics]{par}} 
#' @param imgxSize (integer) width of image when saving to file, see also \code{\link[graphics]{par}} 
#' @param useCol (character or integer) custom colors, see also \code{\link[graphics]{par}} 
#' @param useBord (character) custom histogram elements border color, see also \code{\link[graphics]{par}} 
#' @param isLog (logical) for lin scale signal intensity values where repesentation needs log, assume log2 if \code{TRUE}
#' @param cexSubTi (numerical) subtitle size (expansion factor cex), see also \code{\link[graphics]{par}} 
#' @param cropHist (logical) -not implemented yet- designed for cutting off bars with very low ('insignificant') values
#' @param parDefault (logical) to automatic adjusting par(marg=,cex.axis=0.8), see also \code{\link[graphics]{par}}
#' @param silent (logical) suppress messages
#' @param debug (logical) additonal messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function produces a histogram type graphic (to the ccurrent graphical device)
#' @seealso \code{\link[graphics]{hist}}
#' @examples
#' set.seed(2016); dat1 <- round(c(rnorm(200,6,0.5),rlnorm(300,2,0.5),rnorm(100,17)),2)
#' dat1 <- dat1[which(dat1 <50 & dat1 > 0.2)]
#' histW(dat1,br="FD",isLog=FALSE)
#' histW(log2(dat1),br="FD",isLog=TRUE)
#' 
#' ## quick overview of distributions  
#' layout(partitionPlot(4))
#' for(i in 1:4) histW(iris[,i],isLog=FALSE,tit=colnames(iris)[i])
#' @export
histW <- function(dat, fileName="histW", output="screen", nBars=8, breaks=NULL, tit=NULL, subTi=NULL,xLab=NULL,yLab=NULL,las=NULL,xcex=0.7,
  imgxSize=900, useCol=NULL, useBord=NULL, isLog=TRUE,cexSubTi=NULL,cropHist=TRUE,parDefault=TRUE,silent=FALSE,debug=FALSE,callFrom=NULL) {
  ## make (small) histogram, eg as overview of bandwidth of values in 'dat'
  ## 'isLog'  .. for lin scale SI values where repesentation needs log
  ## 'breaks' .. for (partial) compatibility with hist() : use only for number of breaks (or 'FD'), gets priority over 'nBars'
  ## 'output' .. options for saving image in various formats (set to 'jpg','png',...)
  ## 'subTi'  .. may be FALSE for NOT displaying, or any text, otherwise range
  ## png & co : image height chosen automatically based on number of cols/bars (lower with many bars)
  ## 'cropHist' ..[not implemented yet] designed for cutting off bars with very low ('insignificant') values
  ## 'parDefault' .. to automatic adjusting par(marg=,cex.axis=0.8)
  argN <- deparse(substitute(dat))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="histW")
  opar <- graphics::par(no.readonly = TRUE) 
  on.exit(graphics::par(opar$mar))           # restore only parameters susceptible of changes lateron (to allow successive plotting after layout() )
  on.exit(graphics::par(opar$font.sub))      # restore only parameters susceptible of changes lateron (to allow successive plotting after layout() )
  maxNoXLabels <- NA           # max number of border-values to display
  msg1 <- " NOT sufficent rights to write : "
  if(debug) silent <- FALSE
  if(length(output) <1) output <- "screen" else if(length(output) >1) output <- output[1]
  if(any(output %in% c("png","tiff","tif","jpg","jpeg"))) if(file.access(fileName,mode=0) ==0) {
    if(file.access(fileName, mode=2) <0) stop("File exists ",msg1,fileName) }
  if(length(dat) <1) stop( " 'dat must be vector of numeric values")
  if(is.null(tit)) tit <- paste("Signal Histogram of",argN)
  if(is.null(xLab)) xLab <- "" 
  if(is.null(yLab)) yLab <- paste("Frequency")
  if(is.null(useCol)) useCol <- grDevices::grey(0.8)
  if(is.null(useBord)) useBord <- grDevices::grey(0.7)
  if(length(maxNoXLabels) !=1) maxNoXLabels <- NA
  if(length(breaks)==1) nBars <- breaks
  hist.ini <- graphics::hist(dat, plot=FALSE, breaks=nBars)       
  if(is.character(nBars)) nBars <- length(hist.ini$breaks)
  if(sum(hist.ini$density,na.rm=TRUE) <1) hist.ini$density <- hist.ini$density/sum(hist.ini$density,na.rm=TRUE)
  relCo <- rbind(rel=hist.ini$density, interv=1:length(hist.ini$counts))
  inte <- hist.ini$density
  maxCla <- which(inte==max(inte))
  if(length(inte < nBars) & (sum(inte[maxCla:min(length(inte),maxCla+1)]) > 0.8*sum(inte) |
    sum(inte[maxCla:max(1,maxCla-1)]) > 0.8*sum(inte))) {
    if(!silent) message(fxNa," too much information in 2 neighbour classes (out of ",length(inte),"), trying to increasing number of classes")
    hist.ini <- graphics::hist(dat, plot=FALSE, breaks=nBars+1)              
    relCo <- rbind(rel=hist.ini$density,interv=1:length(hist.ini$counts)) }
  mainHist <- hist.ini
  if(cropHist) {  
     ## crop borders to avoid showing very low bars not yet finished/implemented
    mainHist <- mainHist }
  imgSize <- round(c(imgxSize, imgxSize*(sort(c(0.8, 1.1 + log10(length(mainHist$breaks))/3,1.4))[2])))             ## improve ?
  pl2fi <- FALSE                   # info if plot to file
  tmp <- NULL
  if(output == "png") {tmp <- try(grDevices::png(file=fileName, width=imgSize[1], height=imgSize[2],res=140), silent=TRUE); pl2fi <- TRUE}
  if(output %in% c("jpg","jpeg")) {tmp <- try(grDevices::jpeg(file=fileName, width=imgSize[1], height=imgSize[2],res=140), silent=TRUE); pl2fi <- TRUE}
  if(output %in% c("tif","tiff")) {tmp <- try(grDevices::tiff(file=fileName, width=imgSize[1], height=imgSize[2],res=140), silent=TRUE); pl2fi <- TRUE}
  ## create plot
  if(inherits(tmp, "try-error")) {
    warning(fxNa,"Cannot create file '",output,"' (or other selected format) !")
  } else  {
    nBars <- length(mainHist$breaks)
    ## try to adopt number of x-axis labels displayed to space available
    ch1 <- if(pl2fi) imgSize[1]/150 else graphics::par("pin")[1]            # bring px-size & screen-width to similar measure ..
    nCharLab <- min(max(nchar(c(utils::head(mainHist$breaks), utils::tail(mainHist$breaks)))),5)   # max no of character used in legend (but no more than 5)
    ## how many axis labels can one reasonably display ?
    if(debug) message("Initial ch1 ",signif(ch1,3),"   nCharLab ",nCharLab,"   nBars ",nBars)
    ch1 <- max(2+ ch1 *5 / (1 +nCharLab*xcex), 3)
    if(length(maxNoXLabels)==1 & is.finite(maxNoXLabels)) if(nBars >maxNoXLabels) ch1 <- max(nBars %/% maxNoXLabels, ch1)
    ch2 <- round(nBars/ch1)
    ch1 <- round(ch1) 
    if(debug) message(" final ch1 ",signif(ch1,3),"   nBars ",nBars,"   ch2 ",ch2)
    showBorder <- if(nBars >3 & nBars >ch1) ch2*(1:(nBars/ch2)) else 1:nBars 
    ## define las                                           # x-axis numbers text orientation
    txtLas <- if(length(las) !=1) 1+ (nCharLab >4) else las
      
    if(parDefault) graphics::par(mar=c(0.8 +txtLas, 5, 2.6, 0.8), cex.axis=xcex)                    # mar(bot,le,top,ri)
    cexSubTi <- if(!is.numeric(cexSubTi) | length(cexSubTi) !=1) graphics::par("font.sub")*0.7
    txt2 <- c("Initial values range from  ","  to  ")
    graphics::plot(mainHist, xlab=xLab, ylab=yLab, col=useCol, border=useBord, main=tit, xaxt="n", las=1)
    graphics::segments(hist.ini$breaks[showBorder], 0,hist.ini$breaks[showBorder],  -1*signif(max(hist.ini$counts,na.rm=TRUE)/110,3) ) # x-axis ticks
    ## look for alternatve of making plot, rather as contour-lines for overlay ?
    if(isLog) {
      tx <- signif(c(2^min(dat), 2^hist.ini$breaks[-1*c(1,length(hist.ini$breaks))], 2^max(dat)),2)[showBorder]
      if(length(tx) > length(unique(tx))) tx <- signif(c(2^min(dat), 2^hist.ini$breaks[-1*c(1,length(hist.ini$breaks))], 2^max(dat)),3)[showBorder]
      graphics::mtext(tx, at=hist.ini$breaks[showBorder], col=1, side=1, li=-0.2,cex=0.75, xlab=xLab, las=txtLas)
      txt2 <- paste(txt2[1], signif(2^min(dat,na.rm=TRUE),4), txt2[2], signif(2^max(dat,na.rm=TRUE),5))      # 'initial values range form ..to..'
    } else {
      tx <- signif(hist.ini$breaks[-1*length(hist.ini$breaks)],2)[showBorder]
      if(length(tx) > length(unique(tx))) tx <- signif(c(min(dat), hist.ini$breaks[-1*c(1,length(hist.ini$breaks))], max(dat)),3)[showBorder]
      graphics::mtext(tx, at=hist.ini$breaks[showBorder], col=1, side=1, li=-0.2, cex=0.75,xlab=xLab,las=1)
      txt2 <- paste(txt2[1], signif(min(dat,na.rm=TRUE),4), txt2[2], signif(max(dat,na.rm=TRUE),5)) }
    if(!identical(subTi,FALSE)) graphics::mtext(if(is.null(subTi)) txt2 else subTi,li=-0.4,cex=cexSubTi)
  }
  if(output %in% c("png","tif","tiff","jpg","jpeg")) grDevices::dev.off() }
    
