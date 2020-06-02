#' Add bagplot to existing plot
#' 
#' This function adds a bagplot on an existing (scatter-)plot allowing to highlight the central area of the data.
#' Briefly, a bagplot is a bivariate boxplot, see \href{https://en.wikipedia.org/wiki/Bagplot}{Bagplot}, following the basic idea of a boxplot in two dimensions.
#' Of course, multimodal distributions - if not separated first - may likely lead to mis-interpretation, similarly as it is known for interpreting boxplots.
#' If a group of data consists only of 2 data-points, they will be conected using a straight line.
#' It is recommended using transparent colors to highlight the core part of a group (if only 2 points are available, they will be conected using a straight line),
#' in addition, one could use the option to re-plot all (non-outlyer) points (arguments \code{reCol}, \code{rePch} and \code{reCex} must be used).
#' 
#' @param x (matrix, list or data.frame) main numeric input of data/points to plot
#' @param lev1 (numeric) min content of data for central area (default 0.5 for 50 percent)
#' @param outCoef (numeric) parameter for defining outliers (equivalent to \code{range} in \code{\link[graphics]{boxplot}})
#' @param bagCol (character or integer) color for filling center part of bagplot (default light transparent grey); Note: It is highly suggested to use transparency, otherwise points underneith will be covered
#' @param bagCont (character) color for inner and outer contours of bagplot
#' @param bagLwd (numeric) line width for outer contour, set to \code{NULL} for not diaplaying outer contour (see also \code{\link[graphics]{par}})
#' @param nCore (integer) decide when center should be determined by median or mean: if number of points reach \code{nCore} the median will be used 
#' @param outlCol (character or integer) color for highlighting outlyers (for text and replottig outlyers points), set to \code{NULL} for not highlighting outlyers at all
#' @param outlPch (integer) symbol replotting highlighted outlyers (for text and replottig outlyers points), set to \code{NULL} for not replotting outlyer-points (see also \code{\link[graphics]{par}})
#' @param outlCex (numeric) cex type expansion factor for labels of highlighted outlyers, set to \code{NULL} for not printing (row)names of outlyers (see also \code{\link[graphics]{par}})
#' @param reCol (character or integer) color for replotting (non-outlyer) points, default set to \code{NULL} for not replotting
#' @param rePch (integer) symbol for replotting (non-outlyer) points, default set to \code{NULL} for not re-plotting (see also \code{\link[graphics]{par}})
#' @param reCex (numeric) cex type expansion factor for lfor replotting (non-outlyer) points, default set to \code{NULL} for not replotting
#' @param ctrPch (integer) symbol for shwing group center (see also \code{\link[graphics]{par}})
#' @param ctrCex (numeric) cex type expansion factor for size of group center (see also \code{\link[graphics]{par}})
#' @param ctrCol (character or integer) color for group center symbol
#' @param returnOutL (logical) decide if rownames of (potential) outlyer values should be returned when running the function
#' @param callFrom (character) allow easier tracking of messages produced
#' @param silent (logical) suppress messages
#' @return plot, optional return of matrix with outlyers
#' @seealso  \code{\link{plotPCAw}},  \code{\link[stats]{princomp}}
#' @examples
#' set.seed(2020); dat1 <- matrix(round(rnorm(2000),3),ncol=2); rownames(dat1) <- 1:nrow(dat1)
#' dat1 <- dat1+ 5*matrix(rep(c(0,1,1,0,0,0,1,1),nrow(dat1)/4),byrow=TRUE,ncol=2)
#' col1 <- rgb(red=c(153,90,203,255 ),green=c(143,195,211,125),blue=c(204,186,78,115),
#'   alpha=90,maxColorValue=255)
#' ## suppose we know the repartition into 4 subgroups which we would like to highlight them
#' grp1 <- rep(1:4,nrow(dat1)/4)
#' plot(dat1,col=grey(0.8),xlab="x",ylab="y",las=1,pch=grp1)
#' for(i in 1:4) addBagPlot(dat1[which(grp1==i),],bagCol=col1[i])
#' ## slightly improved
#' library(wrMisc)
#' col2 <- convColorToTransp(col1,255)
#' plot(dat1,col=grey(0.8),xlab="x",ylab="y",las=1,pch=grp1)
#' for(i in 1:4) addBagPlot(dat1[which(grp1==i),],bagCol=col1[i],outlPch=i,
#'   outlCol=col2[i],bagLwd=3)
#' @export
addBagPlot <- function(x,lev1=0.5,outCoef=2,bagCol=NULL,bagCont=bagCol,bagLwd=1.5,nCore=4,outlCol=2,outlPch=NULL,outlCex=0.6,reCol=NULL,rePch=NULL,reCex=NULL,
  ctrPch=NULL,ctrCol=NULL,ctrCex=NULL,returnOutL=FALSE,callFrom=NULL,silent=TRUE) {            # colOutL  colCont=NULL,colOutlP=2,colOutlT=2,
  ##  'x' should be matrix or dataframe (use 1st & 2nd column, ie x & y coord for points) to draw simple bag-plot
  ## 'lev1' gives the min % of points to be included to core (shaded using 'bagCol'), as long as >nCore data-points available
  ## "outliers" are determined similar to boxplots using the 'outCoef'-parameter and then shown in color 'colOutL' and their names may be exported
  ## optional: overall contour (wo outliers) if 'colCont' (=color for contour) given, show center (median) if 'ctrPch' given
  fxNa <- wrMisc::.composeCallName(callFrom,newNa=".addBagPlot")
  msg <- " 'x' must be matrix or data.frame (with at least 1 row and 2 columns)"
  if(is.null(x) <0) stop(msg)
  if(!"matrix" %in% class(x)) x <- try(as.matrix(x))
  if(any(class(x) == "try-error")) stop(msg," - could not transform 'x' using as.matrix()")
  if(length(dim(x)) <2) {
    x <- matrix(x,ncol=length(x))
    if(!silent) message(fxNa," : trying to transform .. since ",msg)
  } else if(length(dim(x)) >2) {x <- as.matrix(x[,,1])}
  if(ncol(x) <2) stop("'x' should have at least 2 columns")
  if(ncol(x) >2 & !silent) {message(fxNa," : only the first 2 columns of data will be used"); x <- x[,1:2]}
  if(length(bagCol) <1) bagCol <- grDevices::rgb(0.1,0.1,0.1,0.1)
  if(length(x) >0) {
    if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
    ctr <- if(nrow(x) < nCore) colMeans(x,na.rm=TRUE) else apply(x,2,stats::median,na.rm=TRUE)                      # overall center : medain if >5 elements
    di <- sqrt((x[,1]-ctr[1])^2 + (x[,2]-ctr[2])^2)                # Euclidean distance to center
    keepX <- if(nrow(x) > 2) di <= stats::median(di,na.rm=TRUE)+ outCoef*diff(stats::quantile(di,c(0.25,0.75))) else rep(TRUE,nrow(x))
    if(sum(keepX) <1) { keepX <- rep(TRUE,nrow(x))
      if(!silent) message("problem defining non-outlyer part of data, keep all")}
    ## define outlyers
    outL <- matrix(x[which(!keepX),],ncol=2)
    if(nrow(outL) >0) rownames(outL) <- rownames(x)[which(!keepX)]
    offS <- if(nrow(x) >1) apply(x,2,function(z) max(abs(range(z,finite=TRUE)),na.rm=TRUE))/70 else x/70               # 
    if(!silent)  {if(nrow(x) > 2) message(fxNa," keep ",sum(keepX)," out of ",nrow(x)," and consider ",
      sum(!keepX)," as outliers") else message(fxNa," too few data, use all columns")}
    if(sum(keepX) < nrow(x)) { x <- x[which(keepX),]
      di <- di[which(keepX)] }
    ## chull around core data
    xCore <- x
    liBag <- di <= stats::quantile(di,lev1,na.rm=TRUE) +min(di,na.rm=TRUE)/100
    if(sum(liBag) < length(liBag)/2.7) liBag <- di <= stats::quantile(di,lev1,na.rm=TRUE) +mean(di,na.rm=TRUE)/10
    if(sum(liBag) >1 & sum(liBag) < nrow(x)) xCore <- x[which(liBag),]
    htps <- grDevices::chull(xCore)
    if(nrow(x) > 2) {
      ## shade core ...    
      graphics::polygon(xCore[c(htps,htps[1]),],col=bagCol,border=bagCont)   
      ## draw outer contour :
      htps2 <- grDevices::chull(x)
      if(length(bagLwd) >0 & !all(is.na(bagCont))) { y <- x[htps2,]; y <- cbind(y,y[c(2:nrow(y),1),])
        graphics::segments(y[,1],y[,2],y[,3],y[,4],col=bagCont,lwd=bagLwd) }
      ## optional replotting of non-outlyer-points
      if(length(reCol) >0 & length(rePch) >0 & length(reCex) >0) graphics::points(x,pch=rePch,col=reCol,cex=reCex)      
    } else if(nrow(x)==2) graphics::lines(x[,1],x[,2],lwd=5,col=bagCol)      # can only connect 2 remaining points by fat line
    ## show group-center (only if ctrPch defined)
    if(is.null(ctrCol)) ctrCol <- bagCol
    if(!is.null(ctrPch)) {            
      if(is.null(ctrCex)) ctrCex <- 1.5
      graphics::points(ctr[1],ctr[2],pch=ctrPch,col=ctrCol,cex=ctrCex) }
    ## highlight outliers
    if(length(outlCol) >0 & nrow(outL) >0) {
      if(length(outlPch) >0) graphics::points(outL,pch=outlPch,col=outlCol)
      if(length(outlCex) >0) graphics::mtext(paste("names of ",sum(!sapply(outL,is.null), na.rm=TRUE),
        " elements looking like potential outlyers were displayed"), cex=0.55, line=-0.8, col=grDevices::grey(0.4))
      if(length(outlCex) >0) graphics::text(outL[,1]+offS[1],outL[,2]+offS[2],col=outlCol, adj=0,cex=outlCex,labels=substr(rownames(outL),1,21))
      }
  }
  if(returnOutL) return( if(is.null(rownames(outL))) which(!keepX) else rownames(outL))
}
  
