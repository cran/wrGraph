#' Find best place on plot for placing legend
#'
#' This function tries to find the best location for placing a legend of a bivariate plot, ie scatter-plot.
#' All 4 corners of the data to plot are inspected for the least occupation by data plotted while displaying the content of \code{sampleGrp}.
#' Alternatively, by setting the argument \code{showLegend} the user-defined legend will be returned 
#'	 
#' @param matr (matrix, list or data.frame) main data of plot
#' @param sampleGrp (character or factor) with this option the text to be displayed in the legend may be taken into consideration for its length
#' @param showLegend (logical or character) decide if \code{matr} should be checked for best location; if \code{showLegend} contains any of the standard legend-location designations (eg 'topleft') it will be used in the output
#' @param suplSpace (numeric) allows to consider extra room taken in legend by symbol and surrounding space, interpreted as n additional characters
#' @param testCorner (integer) which corners should be considered (1=left-top, 2=right-top, right-bottom, left-bottom) 
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @return list with $showL indicating if legend is desired and $loc for the proposition of the best location, $nConflicts gives the counts of conflicts
#' @seealso \code{\link[graphics]{legend}}
#' @examples
#' dat1 <- matrix(c(1:5,1,1:5,5), ncol=2)
#' grp <- c("abc","efghijk")
#' (legLoc <- checkForLegLoc(dat1, grp)) 
#' plot(dat1, cex=3)
#' legend(legLoc$loc, legend=grp, text.col=2:3, pch=1, cex=0.8)
#' @export
checkForLegLoc <- function(matr,sampleGrp=NULL,showLegend=TRUE,suplSpace=4,testCorner=1:4,silent=TRUE,callFrom=NULL){
  ## check 'showLegend' if specified location or logical term (if TRUE, then estimate best location using .bestLegendLoc())
  ## 'matr' .. matrix of data for display, optimal legend location will get determined for these data (if not specifically given)
  ## 'sampleGrp' .. legend-names to be displayed
  ## 'testCorner' specifies places(corners) to be tested: c("topleft","topright","bottomright","bottomleft")
  ## return list with $showL as logical depending if legend should be drawn, and  $loc as location
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="checkFForLegLoc")
  txt <- "'matr' must be matrix or data.frame with at least 1 row and at least 2 columns"
  if(length(dim(matr)) <2) stop(txt)
  if(any(dim(matr)[1:2] < 1:2)) stop(txt)
  .longTxt <- function(x,nLim=8) {    # return longest or almost longest of text-entries in 'x': if more than 'nLim' text-entries look for 87%quantile
    out <- NULL                # initialize
    if(any(nchar(x) >0)) {
      y <- graphics::strwidth(as.character(x), units="figure")
      if(length(x) >nLim) {
        out <- x[order(y, decreasing=TRUE)][round(0.13*length(x))]
      } else out <- x[which.max(y)] } 
  out }
  displLeg <- list(showL=FALSE, loc=NA)
  if(length(showLegend) <1) showLegend <- FALSE
  chLeg <- try(is.logical(showLegend), silent=TRUE)
  if("try-error" %in% class(chLeg)) {            # non-logical entry, check if precise location
    chLoc <- showLegend[1] %in% c("topleft","topright", "top","bottomright", "bottomleft",
      "bottom", "left","right","center")
    if(silent) message(fxNa," argument 'showLegend' was called as TRUE but with  '",showLegend,"', (",if(chLoc) "valid" else "INVALID"," location)")  
    displLeg <- list(showL=chLoc, loc=if(chLoc) showLegend[1] else NA, nConflicts=NA)
  } else {                                                                  # is logical
    if(chLeg) {                                                             # chLeg is TRUE
      sampleGrp <- as.character(wrMisc::naOmit(sampleGrp))
      txLe <- if(length(sampleGrp) >0) .longTxt(sampleGrp) else "sample"
      txLe <- graphics::strwidth(txLe,units="figure") +suplSpace*graphics::strwidth("z", units="figure") 
      nLi <- if(length(sampleGrp) >0) length(unique(sampleGrp)) else 4      # presume 4 groups (if no names given)
      txHi <- (nLi+2.2)*graphics::strheight("1",units="figure")             # need to adjust for extra space towards edge/box of legend
      locCount <- .bestLegendLoc(matr[,1:2],txtLen=txLe,txtHi=txHi,silent=silent,callFrom=callFrom) 
      ## prefer legend on left, add small penalty to right locations to favour lactions at left side (in case of egality)
      if(length(displLeg) >3) locCount[2:3] <- locCount[2:3] +0.1
      best <- which.min(locCount[testCorner])
      if(length(displLeg) >3) locCount[2:3] <- locCount[2:3] -0.1       # re-correct to real counts
      displLeg <- list(showL=TRUE, loc=names(locCount)[best], nConflicts=locCount)
    } else displLeg <- list(showL=FALSE, loc=NA, nConflicts=NA)
  }
  displLeg }
  
#' @export
.bestLegendLoc <- function(dat,txtLen=NULL,txtHi=NULL,silent=TRUE,callFrom=NULL) {
  ## try to find best corner for legend ; 'dat': matrix or df (use 1st & 2nd column, ie x & y coord for points)
  ## 'txtLen' .. text width from graphics::strwidth()
  ## 'txtHi' .. text height from graphics::strheight() (including inter-line)
  if(!(is.matrix(dat) | is.data.frame(dat))) stop(" 'dat' must be matrix or data.frame (with at least 2 columns)")
  if(ncol(dat) <2) stop("'dat' should have at least 2 columns")
  if(!is.numeric(txtHi)) stop("argument 'txtHi' should be numeric (of length 1) !")
  if(length(txtLen) >1 | length(txtHi) >1) {txtLen <- txtLen[1]; txtHi <- txtHi[1]; warning("truncating 'txtLen' and 'txtHi' to length 1 !")}
  if(is.na(txtLen) | is.na(txtHi) |  min(c(txtLen,txtHi)) < 0)  stop("arguments 'txtLen' and 'txtHi' should not be NA or negative !")
  raX <- range(dat[,1], na.rm=TRUE)
  raY <- range(dat[,2], na.rm=TRUE)
  raX <- c(raX, abs(raX[2] -raX[1]))
  raY <- c(raY, abs(raY[2] -raY[1]))
  ocX <- txtLen/diff(graphics::par("plt")[1:2])       # fraction of text occupied in x
  ocY <- txtHi /diff(graphics::par("plt")[3:4])       # fraction of text occupied in y
  locCount <- c(
    topleft    =sum(dat[,1] < raX[1] + ocX*raX[3] & dat[,2] > raY[2] - ocY*raY[3], na.rm=TRUE) ,
    topright   =sum(dat[,1] > raX[2] - ocX*raX[3] & dat[,2] > raY[2] - ocY*raY[3], na.rm=TRUE) ,
    bottomright=sum(dat[,1] > raX[2] - ocX*raX[3] & dat[,2] < raY[1] + ocY*raY[3], na.rm=TRUE) ,
    bottomleft =sum(dat[,1] < raX[1] + ocX*raX[3] & dat[,2] < raY[1] + ocY*raY[3], na.rm=TRUE) )
  if(!silent) { message(callFrom,"  txtLen ",signif(txtLen,3))
    graphics::abline(v=c(raX[1] + ocX*raX[3], raX[2] - ocX*raX[3]), lty=2)
    graphics::abline(h=c(raY[2] - ocY*raY[3], raY[1] + ocY*raY[3]), lty=2)
    message("  locCount: ", paste(locCount,collapse="  ")) }
  locCount }  
     
