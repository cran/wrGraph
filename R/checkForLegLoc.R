#' Find best place on plot for placing legend
#'
#' This function tries to find the best location for placing a legend of a bivariate plot, ie scatter-plot.
#' All 4 corners of the data to plot are inspected for the least occupation by data plotted.
#'	 
#' @param matr (matrix, list or data.frame) main data of plot
#' @param showLegend (logical or character) decide if \code{matr} should be checked for best location; if \code{showLegend} contains any of the standard legend-location designations (eg 'topleft') it will be used in the output
#' @param sampleGrp (character or factor) with this option the text to be displayed in the legend may be taken into consideration for its length
#' @return list with $showL indicating if legend is desired and $loc for the proposition of the best location
#' @seealso \code{\link[graphics]{legend}}
#' @examples
#' dat1 <- matrix(c(1:5,1,1:5,5),ncol=2)
#' grp = c("abc","efghijk")
#' (legLoc <- checkForLegLoc(dat1,grp)) 
#' legend(legLoc$loc,legend=grp,text.col=2:3,pch=1,cex=0.8)
#' @export
checkForLegLoc <- function(matr,sampleGrp=NULL,showLegend=TRUE){
  ## check 'showLegend' if specified location or logical term (if TRUE, then estimate best location using .bestLegendLoc())
  ## 'matr' .. matrix of data for display, optimal legend location will get determined for these data (if not specifically given)
  ## 'sampleGrp' .. legend-names to be displayed
  ## return list with $showL as logical depending if legend should be drawn, and  $loc as location
  displLeg <- list(showL=FALSE,loc=NA)
  if(is.character(showLegend)) {if(any(c("T","F","TRUE","FALSE") %in% showLegend)) displLeg[[1]] <- as.logical(showLegend[1]) else {
    if(showLegend %in% c("topleft","topright", "top","bottomright", "bottomleft",
      "bottom", "left","right","center")) displLeg <- list(showL=TRUE,loc=showLegend[1])}} else {
    if(is.logical(showLegend)) displLeg[[1]] <- showLegend[1]}
  if(displLeg[[1]]) { if(is.na(displLeg[[2]])) displLeg[[2]] <- .bestLegendLoc(
      matr[,1:2],txtLen=min(0.25,round(0.06+max(nchar(as.character(sampleGrp)))/110,2)),
      legLen=min(0.7,round(0.04+length(levels(sampleGrp))/70,2)))}
  displLeg }

#' @export
.bestLegendLoc <- function(dat,testCorner=1:4,txtLen=0.15,legLen=0.3) {
  ## try to find best corner for legend ; 'dat': matrix or df (use 1st & 2nd column, ie x & y coord for points)
  ## 'testCorner' specifies places(corners) to be tested: c("topleft","topright","bottomright","bottomleft")
  if(!(is.matrix(dat) | is.data.frame(dat))) stop(" 'dat' must be matrix or data.frame (with at least 2 columns)")
  if(ncol(dat) <2) stop("'x' should have at least 2 columns")
  if(!is.numeric(txtLen) | !is.numeric(legLen)) stop("arguments 'txtLen' and 'legLen' should be numeric (of length 1) !")
  if(length(txtLen) >1 | length(legLen) >1) {txtLen <- txtLen[1]; legLen <- legLen[1]; warning("truncating 'txtLen' and 'legLen' to length 1 !")}
  if(is.na(txtLen) | is.na(legLen) | min(c(txtLen,legLen)) <= 0 | max(c(txtLen,legLen)) >= 1)  stop("arguments 'txtLen' and 'legLen' should have values >0 and <1 !")
  raX <- range(dat[,1],na.rm=TRUE)
  raY <- range(dat[,2],na.rm=TRUE)
  locCount= as.matrix(data.frame(
    topleft    =sum(dat[,1] < raX[1 ]+ txtLen*abs(diff(raX)) & dat[,2] > raY[2] - legLen*abs(diff(raY)),na.rm=TRUE) ,
    topright   =sum(dat[,1] > raX[2] - txtLen*abs(diff(raX)) & dat[,2] > raY[2] - legLen*abs(diff(raY)),na.rm=TRUE) ,
    bottomright=sum(dat[,1] > raX[2] - txtLen*abs(diff(raX)) & dat[,2] < raY[1] + legLen*abs(diff(raY)),na.rm=TRUE) ,
    bottomleft =sum(dat[,1] < raX[1] + txtLen*abs(diff(raX)) & dat[,2] < raY[1] + legLen*abs(diff(raY)),na.rm=TRUE)  ))
  locCount <- locCount[,testCorner]
  names(locCount)[which(locCount==min(locCount))[1]]  }
     
