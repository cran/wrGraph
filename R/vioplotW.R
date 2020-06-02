#' Violin-plots version W
#'
#' This function allows generating violin-plots using a variety of input formats and offers additional options for colors.   
#' Code from \href{https://CRAN.R-project.org/package=vioplot}{vioplot} was modified based on a omment on
#' \href{http://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html}{Nabble} and 
#' further extended. The package \href{https://CRAN.R-project.org/package=vioplot}{vioplot} does not have to be installed for using this function. 
#' Main input data may now be matrix of list of multiple data-elements (entries may be of variable length),
#'	individual colors for different sets of data or color-gradients can be specified, and display n per set of data was integtated. 
#' Note : Arguments have to be given with full names, lazy evaluation of arguments will not work properly with this function (since '...' is used to capture additional data-sets).    
#' @param x (matrix, list or data.frame) data to plot, or first series of data
#' @param ... (numeric) additional sets of data to plot 
#' @param finiteOnly (logical) eliminate non-finite elements to avoid potential errors (eg when encountering \code{NA})
#' @param range (numerical) custom range  
#' @param hh (numeric) a vector of length one, two or three, defining the smoothing parameter (standard deviation to kernel function, if omited anormal optimal smoothing parameter is used); equivalent to argument \code{h} in \code{\link[vioplot]{vioplot}}; see also \code{\link[sm]{sm.density}} 
#' @param ylim (numeric, length=2) custom limit on y-axis, see also \code{\link[graphics]{par}} 
#' @param nameSer (character) custom label
#' @param horizontal (logical) orientation of plot
#' @param col (character or integer) for gradients may be custom colors or  'rainbow', 'grayscale', 'Spectral' or 'Paired' 
#' @param border (character) custom color for figure border
#' @param lty (integer) type of line(s) (see also \code{\link[graphics]{par}}) 
#' @param las (integer) orientation of axis labels (see also \code{\link[graphics]{par}}) 
#' @param tit (character) custom title to figure
#' @param lwd (integer) width of line(s) (see also \code{\link[graphics]{par}}) 
#' @param rectCol (character) color of rectangle
#' @param colMed (character or integer) color for median (or average) point to show center of group
#' @param pchMed (integer) symbol type for center  (see also \code{pch} in \code{\link[graphics]{par}}) 
#' @param at (numeric) custom locoation of data-series names
#' @param add (logical) add to existing plot if \code{TRUE}
#' @param wex (integer) relative expansion of the violin
#' @param drawRect (logical) if \code{TRUE}, draw box
#' @return figure only
#' @seealso (the original :) \code{\link[vioplot]{vioplot}}
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1),ncol=6,dimnames=list(paste("li",1:50,sep=""),
#'   letters[19:24]))
#' vioplotW(datT6)
#' ## variable number of elements (each n is displayed)
#' datT6b <- apply(datT6,2,function(x) x[which(x <5)])
#' vioplotW(datT6b,col="Spectral")
#' @export
vioplotW <- function(x, ..., finiteOnly=TRUE,range=1.5, hh=NULL, ylim=NULL, nameSer=NULL, horizontal=FALSE, 
	col="rainbow", border="black", lty=1,tit=NULL,las=1, lwd=1, rectCol="black", colMed="white", pchMed=19,
    at=0, add=FALSE, wex=1, drawRect=TRUE, silent=FALSE, callFrom=NULL) {
  ## wr modification of vioplot() for violin-plots, modified &extended based on http://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html
  ## 'col'.. individual fill-colors, may specify gradients (default 'rainbow','grayscale','Spectral','Paired')
  ## default display (short)column names and n
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="vioplotW")
  fxArg <- c("x","finiteOnly","range","h","ylim","nameSer","horizontal","col","border",
    "lty","tit","lwd","rectCol","colMed","pchMed","at","add","wex","drawRect") 
  argN <- c(x=deparse(substitute(x)), sup=deparse(substitute(...)))
  colNx <- if(length(dim(x)) >1) colnames(x) else argN[1]
  datas <- wrMisc::asSepList(list(x, ...),asNumeric=TRUE,fxArg=fxArg)
  if(finiteOnly) for(i in  1:length(datas)) {
    chFini <-  is.finite(datas[[i]])
    if(any(!chFini)) datas[[i]] <- datas[[i]][which(chFini)] }
  ## prepare input data
  n <- length(datas)
  if(length(col) <n) col <- rep(col,n)[1:n]
  if(missing(at)) at <- 1:n
  upper <- lower <- nByCol <- vector(mode="numeric", length=n)
  q1 <- q3 <- vector(mode="numeric", length=n)
  med <- vector(mode="numeric", length=n)
  base <- vector(mode="list", length=n)
  height <- vector(mode="list", length=n)
  baserange <- c(Inf, -Inf)
  if("rainbow" %in% col) col <- grDevices::rainbow(round(n*1.07))[1:n]
  if("grayscale" %in% col) col <- grDevices::gray.colors(n)
  if("Spectral" %in% col) col <- RColorBrewer::brewer.pal(min(n,11),"Spectral")
  if("Paired" %in% col)  col <- RColorBrewer::brewer.pal(12,"Paired")[c(5:6,1:4,7:12)][1:min(n,12)]
  args <- list(display = "none")
  if(!(is.null(hh))) args <- c(args, hh = hh)
  for(i in 1:n) {
    data <- as.numeric(as.matrix(datas[[i]]))    # faster than unlist()
    dataMin <- min(data,na.rm=TRUE)
    dataMax <- max(data,na.rm=TRUE)
    q1[i] <- stats::quantile(data, 0.25,na.rm=TRUE)
    q3[i] <- stats::quantile(data, 0.75,na.rm=TRUE)
    med[i] <- stats::median(data,na.rm=TRUE)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, dataMax)
    lower[i] <- max(q1[i] - range * iqd, dataMin)
    est.xlim <- c(min(lower[i], dataMin), max(upper[i], dataMax))
    smDensity <- sm::sm.density 
    smout <- do.call("smDensity", c(list(data, xlim=est.xlim), args))       
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
    nByCol[i] <- length(data)
  }
  if(!add) {
    xlim <- if(n == 1) at + c(-0.5, 0.5) else range(at) + min(diff(at))/2 * c(-1, 1)
    if(is.null(ylim)) ylim <- baserange
  }
  label <- if(is.null(nameSer)) names(datas) else nameSer
  if(is.null(label) | all(is.na(label))) label <- 1:n
  boxwidth <- 0.05 * wex
  nDisp <- sapply(datas,function(x) sum(!is.na(x)))
  nDisp <- if(length(unique(nDisp)) >1) paste(rep(c("n=",""), if(n>6) c(1,n-1) else c(n,0)),nDisp,sep="") else c(paste("n=",nDisp[1]),rep("",n-1))
  if(!add) graphics::plot.new()
  if(!horizontal) {
    if(!add) {
        graphics::plot.window(xlim=xlim, ylim=ylim)
        graphics::axis(2,las=las)
        graphics::axis(1, at=at, label=label,las=las)
        if(!is.null(tit)) graphics::title(main=tit)
    }
    graphics::box()
    for(i in 1:n) {
        graphics::polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
            c(base[[i]], rev(base[[i]])), col=col[i], border=border, lty=lty, lwd=lwd)
        if(drawRect) {
            graphics::lines(at[c(i, i)], c(lower[i], upper[i]), lwd=lwd, lty=lty)
            graphics::rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, q3[i], col=rectCol)
            graphics::points(at[i], med[i], pch = pchMed, col = colMed) }
    }
    graphics::mtext(nDisp,at=1:n,side=3,cex=0.7,line=-1)
  } else {
    if(!add) {
        graphics::plot.window(xlim = ylim, ylim = xlim)
        graphics::axis(1)
        graphics::axis(2, at = at, label = label)
    }
    graphics::box()
    for(i in 1:n) {
        graphics::polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
            rev(at[i] + height[[i]])), col=col[i], border=border, lty=lty, lwd=lwd)
        if (drawRect) {
            graphics::lines(c(lower[i], upper[i]), at[c(i, i)], lwd=lwd, lty=lty)
            graphics::rect(q1[i], at[i] -boxwidth/2, q3[i], at[i] + boxwidth/2, col=rectCol)
            graphics::points(med[i], at[i], pch=pchMed, col=colMed) }
    }
    graphics::mtext(nDisp, at=1:n, side=2, cex=0.7, line=-1)
  }
  invisible(list(upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}
    
