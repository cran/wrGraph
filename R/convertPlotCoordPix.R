#' Convert points of plot to coordinates in pixels 
#' 
#' This function allows conversion the plotting positions ('x' and 'y' coordiantes) of points in a given plot into coordiantes in pixels (of the entire plotting region).
#' It was designed to be used as coordinates in an html file for mouse-over interactivity (display of names of points and links).
#' Of course, the size of the plotting region is crucial and may not be changed afterwards (if the plot is not written to file using \code{png} etc).
#' In turn the function \code{\link{mouseOverHtmlFile}} will use the pixel-coordiantes, allowing to annotate given points of a plot for mouse-over interactive html.
#'	 
#' @param x (numeric) initial plotting coordinates on x-axis, names of vector - if available- will be used as IDs
#' @param y (numeric) initial plotting coordinates on y-axis 
#' @param useMar (numeric,length=4) margins defined with plot, see also \code{\link[graphics]{par}} 
#' @param plotDim (integer, length=2) dimension of the plotting device in pixels, see also \code{\link[graphics]{par}}  
#' @param plotRes (integer) resoltion of plotting device, see also \code{\link[graphics]{par}}
#' @param fromTop (logical) toggle if poordinates should start from top 
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @return matrix with x- and y-coordinates in pixels
#' @seealso  \code{\link{mouseOverHtmlFile}}
#' @examples
#' df1 <- data.frame(id=letters[1:10],x=1:10,y=rep(5,10),mou=paste("point",letters[1:10]),
#'   link=file.path(tempdir(),paste(LETTERS[1:10],".html",sep="")),stringsAsFactors=FALSE)  
#' ## Typically one wants to get pixel-coordinates for plots written to file.
#' ## Here we'll use R's tempdir, later you may want to choose other locations
#' pngFile <- file.path(tempdir(),"test01.png")
#' png(pngFile,width=800, height=600,res=72)
#' ## here we'll just plot a set of horiontal points at default parameters ...
#' plot(df1[,2:3],las=1,main="test01")
#' dev.off()
#' ## Note: Special characters should be converted for proper display in html during mouse-over
#' library(wrMisc)
#' df1$mou <- htmlSpecCharConv(df1$mou)
#' ## Let's add the x- and y-coordiates of the points in pixels to the data.frame
#' df1 <- cbind(df1,convertPlotCoordPix(x=df1[,2],y=df1[,3],plotD=c(800,600),plotRes=72))
#' head(df1)
#' ## using mouseOverHtmlFile() one could now make an html document with interactive 
#' ## display of names and clockable links to the coordinates we determined here ...
#' @export
convertPlotCoordPix <- function(x,y,useMar=c(6.2,4,4,2),plotDim=c(1400,800),plotRes=100,fromTop=TRUE,callFrom=NULL,silent=FALSE){
  ## convert point-coordinates 'x' and 'y' (of plot) in pixel coordinates  (eg for use with mouse-over tip in Html)
  ## return matrix with 2 columns :  pxX & pxY .. coresponding pixel coordinates (with length(x) rows)
  ## expecting 'useMar' as margins given in par() (as bottom,left,top,right) in inch
  ## 'plotDim' should be image size for png-device in pixels (width, height), 'plotRes' should be 'plotRes' for png()
  ## (default) 'fromTop'=TRUE : coordinates will be counted from top border of image
  ##  like png(plotDim[1],plotDim[2,plotRes=plotRes]
  fxNa <- wrMisc::.composeCallName(callFrom,newNa=" convertPlotCoordPix")
  msg <- " arguments 'x' & 'y' should be finite numeric and of same length ! "
  if(any(length(x) != length(y),length(x) <1,length(y) <1)) stop(msg)
  isFinit <- cbind(x=is.finite(x),y=is.finite(y))
  if(any(!isFinit)) {
    if(all(!isFinit[,1]) | all(!isFinit[,2])) stop(msg)
    if(any(!isFinit[,1])) x[which(!isFinit[,1])] <- sum(x[which(isFinit[,1])])/length(isFinit[,1])
    if(any(!isFinit[,2])) y[which(!isFinit[,2])] <- sum(y[which(isFinit[,2])])/length(isFinit[,2]) }
  msg0 <- "numeric vector of length"
  msg <- paste("'plotDim' should be",msg0,"2 (x & y dimension in pixel)")
  if(length(plotDim) !=2) stop(msg)
  msg <- paste("'useMar' should be",msg0,"4;  'plotRes'",msg0,"1")
  if(any(length(useMar) !=4,length(plotRes) !=1)) stop(msg)
  if(any(!is.finite(useMar),!is.finite(plotRes))) stop(msg)
  if(is.null(names(x))) {
    if(!is.null(names(y))) { names(x) <- names(y) }}
  msg2 <- " 'useMar' should be vector of 4 values (as bottom,left,top,right)"
  msg3 <- " 'plotDim' should be vector of 2 values (as width & hight of png)"
  if(!all(is.finite(useMar)) | length(useMar) !=4) stop(msg2)
  if(any(useMar <0)) stop(msg2)
  if(!all(is.finite(plotDim)) | length(plotDim) !=2) stop(msg3)
  if(any(plotDim <= 0)) stop(msg3)
  ##
  if(!silent) message(fxNa,"supl params : mar=",paste(useMar,collapse=","),"  dim=",
    plotDim[1],"x",plotDim[2],"   res=",plotRes)
  out <- .predPointsPix(x,y,dimPng=plotDim,res=plotRes,marg=useMar,fromTop=fromTop)
  if(any(!isFinit[,1])) out[which(!isFinit[,1]),1] <- NA
  if(any(!isFinit[,2])) out[which(!isFinit[,2]),2] <- NA
  out }

#' @export
.determFigMargPix <- function(marg,res){
  ## estimate size/distance of margin to edge of image (png) in pixel (return numeric vector)
  ## in case of 4 margin values it is assumed that these follow the order as in par() ie c(bottom,left,top,right)
  ## marg .. distance of margin in inch (as given in par(mar=c()),  )
  ## res .. resolution of image (png)
  sloC <-  -0.3735 +0.2034*res                            # based on optimization wr12may15
  offSC <- -0.2 -4.714e-12*(res^5)
  dis <- floor(offSC + sloC*(marg))
  dis <- dis - if(length(dis) ==4) c(1,2,1,2) else 1      # additional correction (based on 800x600 at res=100)
  dis[dis <0] <- 0
  if(length(dis) ==4) names(dis) <- c("bottom","left","top","right")
  dis }

#' @export
.predPointsPix <- function(x,y,dimPng,res,marg,fromTop=TRUE,scExt=0.04,displ=FALSE){
  ## predict & return pixel location of points of (default) plot() (which uses 4% extension of scales)
  ## return 2 col matrix with columns 'xPix' and 'yPix' (with length(x) rows)
  ## 'x' & 'y' .. initial coordiantes for plot (to be made separately)
  ## 'pngDim' as width & hight of png
  ## 'res' as resultion of png
  ## 'marg' as  par(mar=c(parMargTop, parMargLe, parMargBot, parMargRi) )
  ## 'fromTop' ... default counting in html is from top
  ## 'scExt' .. extending scale at 4%
  ## note: may be imprecise in case of x or y with all same values 
  plPx <- matrix(abs(.determFigMargPix(marg,res)[c(2,4,3,1)]+c(0,-1*dimPng[1],0,-1*dimPng[2])),ncol=2,dimnames=list(c("start","end"),c("x","y")))
  if(displ) graphics::plot(x,y,las=1)
  plaRa <- abs(c(diff(range(x)),diff(range(y))))
  ##  extRa  .. pos in pix of corners (le,ri,bot,top) of plot in plot-coord
  extRa <- c(range(x)+plaRa[1]*scExt*c(-1,1), range(y)+plaRa[2]*scExt*c(-1,1))    # value of extended x-scal (new min & max), ext y-scale (new min&max)
  if(displ) {graphics::points(extRa[1],extRa[3],pch=3,col=2); graphics::points(extRa[2],extRa[4],pch=3,col=3)}    # 
  if(length(unique(x))==1) extRa[1:2] <- range(pretty(unique(x)*c(0.55,1))) +plaRa[1]*scExt*c(-1,1)       
  if(length(unique(y))==1) extRa[3:4] <- range(pretty(unique(x)*c(0.55,1))) +plaRa[2]*scExt*c(-1,1)  
  xIncr <- diff(plPx[,1])/diff(extRa[1:2])     # how many pix per x-unit
  yIncr <- diff(plPx[,2])/diff(extRa[3:4])     # how many pix per y-unit
  xPix <- plPx[1,1] +(x-extRa[1])*xIncr                           
  yPix <- if(fromTop) plPx[1,2] +diff(plPx[,2])- (y-extRa[3])*yIncr  else {
    dimPng[2] -plPx[2,2] +(y-extRa[3])*yIncr }       # calculate as counted from top or bottom of png-image
  round(cbind(xPix,yPix)) }
   
