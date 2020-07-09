#' Add histogram to existing plot  
#'
#' Add histogram at pleace of legend using colors from 'colorRamp'.
#'
#' @param x (numeric) main input/component of plot
#' @param colRamp (character or integer) set of colors, default is rainbow-like
#' @param location (character) for location of histogram inside existing plot (may be 'br','bl','tl','tr','bottomright', 'bottomleft','topleft','topright') 
#' @param legTit (character, length=1) optional title for histogram-insert
#' @param border (logical) draw border around histogram-insert
#' @param cex (numeric) expansion factor (see also \code{\link[graphics]{par}})
#' @param srt (numeric) angle for histogram text labels (90 will give vertical label) (see also \code{\link[graphics]{par}})
#' @param offS (\code{NULL} or numeric, length=5) fine-tuning of where histogram-insert will be placed and how elements therein are ditributed 
#'  (default c(xOff=0.2,yOff=0.25,leftOffS=0.05, upperBarEnd=1.05,txtOff=0.02),  
#'  1st and 2nd determine proportio of insert relative to entire plotting region, 3rd defines space left on bottom for text, 
#'  4th if bars hit ceiling of insert or proportion to leave, 5th for shifting text towards top when turned other than 90 degrees ) 
#' @param border (logical) decide of draw gray rectangle or not around legend
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return figure
#' @examples
#' dat <- rnorm(90); plot(dat)
#' legendHist(dat, col=1:5)
#' @export
legendHist <- function(x,colRamp=NULL,location="bottomright",legTit=NULL,cex=0.7,srt=67,offS=NULL,border=TRUE,silent=FALSE,callFrom=NULL){
  ## add histogram instead of legend using colors from 'colorRamp', so far as bottomright
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="legendHist")
  if(is.null(colRamp)) {
    colRamp <- cbind(red=c(141,72,90,171, 220,253,244,255),green=c(129,153,194,221, 216,174,109,0),blue=c(194,203,185,164, 83,97,67,0))   # rainbow-like gradient    
    colRamp <- grDevices::rgb(red=colRamp[,1],green=colRamp[,2],blue=colRamp[,3],maxColorValue=255)}
  cutInt <- wrMisc::cutToNgrp(x,colRamp,NAuse=FALSE)
  tmp <- graphics::par("usr")
  .sw <- function(location,x,offSet) {     # function for setting corner-coordinates according of 'x'
    ## uses only 1st & 2nd val of offSet
    if(length(location) !=1) location <- "bottomright"
    location <- tolower(as.character(location))
    if(any(c("bottomri","botri","bori","br") %in% location)) location <- "bottomright"
    if(any(c("bottomle","botle","bole","bl") %in% location)) location <- "bottomleft"
    if(any(c("topri","tori","tr") %in% location)) location <- "topright"
    if(any(c("toptle","tole","tl") %in% location)) location <- "topleft"
    ## make int %value of offSet absolute :
    offSet[1:2] <- offSet[1:2]*c(x[2]-x[1],x[4]-x[3]) 
    switch(location,
      bottomright= c(xMin=x[2]-offSet[1],x[2:3],yMax=x[3]+offSet[2]),
      bottomleft= c(xMin=x[1],xMax=x[1]+offSet[1],x[3],yMax=x[3]+offSet[2]),
      topleft= c(xMin=x[1],xMax=x[1]+offSet[1],x[4]-offSet[2],yMax=x[4]),
      topright= c(xMin=x[2]-offSet[1],x[2],x[4]-offSet[2],yMax=x[4]))} 
  if(length(offS) <5) offS <- c(xOff=0.2,yOff=0.25,leftOffS=0.05,upperBarEnd=if(length(legTit) >0) 1.15 else 1.05, txtOff=0.02)       # %of x-range, % of y-range
  legCor <- .sw(location,tmp,offS)
  legCor <- c(legCor, legCor[3]+offS[2]*(legCor[4]-legCor[3]))     # add 5th component: bottom y axis start for Hist (shifted by offS[2] %
  legWi <- abs(c(legCor[2] -legCor[1],legCor[4] -legCor[3]))       # width of insert on x and y
  br <- seq(legCor[1] +legWi[2]*offS[3],legCor[2] -legWi[2]*0,length.out=length(colRamp)+1)   # define breaks on x-axis, uses offS[3], so far nn offset to right
  chBr <- br > legCor[2]
  if(any(chBr)) br[length(br)] <- legCor[2]
  mids <- br[2]-br[1]                      # $breaks on x
  mids <- br[-length(br)]+mids/2           # midPoints on x
  tb <- table(cutInt$grouped)              #      
  tb <- legCor[5]+ (legCor[4]-legCor[5])*tb/(max(tb,na.rm=TRUE)*offS[4])      # upper end of bars & factor offS[4] typically 1.05
  if(border) graphics::rect(legCor[1],legCor[3],legCor[2],legCor[4],col=grDevices::rgb(1,1,1,0.5),border=grDevices::grey(0.7))  # broder of 'legend'
  srtCo <- c((90 -srt)*(br[2] -br[1])/90, (legCor[4]-legCor[3])*offS[5] )     # supl offset for text-labels, uses factor offS[5]
  for(j in 1:length(colRamp)) {
    graphics::text(mids-srtCo[1],rep(srtCo[2]+legCor[3],length(mids)),paste(">",signif(cutInt$legTxt[,1],3)),adj=0,srt=srt,cex=cex) 
    graphics::rect(br[j],legCor[5],br[j+1],tb[j],col=colRamp[j],border=border)                                               # main histogram
    if(length(legTit) >0) graphics::text(mean(legCor[1:2]),max(tb,na.rm=TRUE)+(legCor[2]-legCor[1])*0.03,legTit[1],cex=cex)  # histogram-insert title
    }
}
   
