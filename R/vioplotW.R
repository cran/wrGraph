#' Violin-plots version W
#'
#' This function allows generating \href{https://en.wikipedia.org/wiki/Violin_plot}{Violin plots}) using a variety of input formats and offers additional options for colors. 
#' Main input may be multiple vectors, a matrix or list of multiple data-elements (entries may be of variable length),
#'	individual colors for different sets of data or color-gradients can be specified, and the display of n per set of data was integtated 
#' (based on an inspiration from  \href{http://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html}{Nabble}). 
#' It is also possible to plot pairwise half-violins for easier pairwise-comparisons (using \code{halfViolin="pairwise"}).
#' Many arguments are kept similar to \href{https://CRAN.R-project.org/package=vioplot}{vioplot} (this package is not required for this function.
#' Note : Arguments have to be given with full names, lazy evaluation of arguments will not work properly with this function (since '...' is used to capture additional data-sets).    
#' Note : \href{https://CRAN.R-project.org/package=vioplot}{vioplot} offers better options for plotting formulas 
#'  
#' @param x (matrix, list or data.frame) data to plot, or first series of data
#' @param ... (numeric) additional sets of data to plot 
#' @param finiteOnly (logical) eliminate non-finite elements to avoid potential errors (eg when encountering \code{NA})
#' @param halfViolin (logical or character) decide with \code{TRUE} or \code{FALSE} if full or only half of violins should be plotted, if "pairwise" always 2 data-sets will be plotted back-to-back
#' @param hh (numeric, length <4) smoothing parameter (standard deviation to kernel function, if omited anormal optimal smoothing parameter is used); equivalent to argument \code{h} in \code{\link[vioplot]{vioplot}}; see also \code{\link[sm]{sm.density}} 
#' @param ylim (\code{NULL} or numeric, length=2) custom limit on y-axis, see also \code{\link[graphics]{par}} 
#' @param nameSer (character) custom label for data-sets or columns (length must match number of data-sets)
#' @param horizontal (logical) orientation of plot
#' @param col (character or integer) custom colors or gradients like 'rainbow', 'grayscale', 'heat.colors', 'topo.colors', 'Spectral' or 'Paired',  or you may use colors made by the package \href{https://CRAN.R-project.org/package=colorRamps}{colorRamps}  
#' @param border (character) custom color for figure border
#' @param las (integer) orientation of axis labels (see also \code{\link[graphics]{par}}) 
#' @param lty (integer) line-type for linear regression line (see also \code{\link[graphics]{par}}) 
#' @param tit (character) custom title to figure
#' @param lwd (integer) width of line(s) (see also \code{\link[graphics]{par}}) 
#' @param rectCol (character) color of rectangle
#' @param at (numeric) custom locoation of data-series names, ie the points at which tick-marks are to be drawn, will be passed to \code{\link[graphics]{axis}}, it's length ust match the number of data-sets
#' @param add (logical) add to existing plot if \code{TRUE}
#' @param wex (integer) relative expansion factor of the violin
#' @return figure only
#' @seealso the package \href{https://CRAN.R-project.org/package=vioplot}{vioplot}, \code{\link[sm]{sm}} is used for the density estimation
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @examples
#' set.seed(2013)
#' dat6 <- matrix(round(rnorm(300)+3,1), ncol=6, 
#' 	 dimnames=list(paste("li",1:50,sep=""), letters[19:24])) 
#' vioplotW(dat6)
#' ## variable number of elements (each n is displayed)
#' dat6b <- apply(dat6,2,function(x) x[which(x <5)])
#' dat6b[[4]] <- dat6b[[4]][dat6b[[4]] <4]
#' vioplotW(dat6b,col="Spectral")
#' vioplotW(dat6b,col="Spectral",halfViolin="pairwise",horizontal=TRUE)
#' vioplotW(dat6b,col="Spectral",halfViolin="pairwise",horizontal=FALSE)
#' @export
vioplotW <- function(x, ..., finiteOnly=TRUE, halfViolin=FALSE, hh=NULL, ylim=NULL, nameSer=NULL, horizontal=FALSE, col="rainbow", 
	border="black", lty=1, tit=NULL, las=1, lwd=1, rectCol="black", at=0, add=FALSE, wex=1, silent=FALSE, callFrom=NULL) {
  ## wr variant for violin-plots, inspired&extended based on http://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html
  ## 'col'.. individual fill-colors, may specify gradients (default 'rainbow','grayscale','Spectral','Paired')
  ## default display (short)column names and n
  doPlot <- TRUE
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="vioplotW")
  rangeVa <- 1.5                                                   # how much to extent fitted distribution beyond real min/max (if
  argN <- c(x=deparse(substitute(x)), sup=deparse(substitute(...)))
  if(length(x) <1) { warning(fxNa,"input 'x' seems empty !"); doPlot <- FALSE} 
  chSm <- try(find.package("sm"), silent=TRUE)
  if("try-error" %in% class(chSm)) { doPlot <- FALSE
    message("Cannot find package 'sm' which is needed for function ",fxNa)} 
  fxArg <- c("x","finiteOnly","hh","ylim","nameSer","horizontal","col","border",
    "lty","tit","lwd","rectCol","colMed","pchMed","halfViolin","at","add","wex","drawRect") 
  colNx <- if(length(dim(x)) >1) colnames(x) else argN[1]
  datas <- wrMisc::asSepList(list(x, ...), asNumeric=TRUE, fxArg=argN,callFrom=fxNa)
  n <- length(datas)
  if(doPlot ){  
    if("rainbow" %in% col) col <- grDevices::rainbow(round(n*1.08))[1:n]
    if("grayscale" %in% col) col <- grDevices::gray.colors(n)
    if("heat.colors" %in% col) col <- grDevices::heat.colors(n)
    if("topo.colors" %in% col) col <- grDevices::topo.colors(n)
    chPaCol <- try(find.package("RColorBrewer"), silent=TRUE)
    if("try-error" %in% class(chPaCol) & any(c("Spectral","Paired") %in% col)) stop("Cannot find package 'RColorBrewer' which is needed for your choice of argment 'col' !")   
    if("Spectral" %in% col & !"try-error" %in% class(chPaCol)) col <- RColorBrewer::brewer.pal(min(n,11),"Spectral")
    if("Paired" %in% col & !"try-error" %in% class(chPaCol))  col <- RColorBrewer::brewer.pal(12,"Paired")[c(5:6,1:4,7:12)][1:min(n,12)]
    if(finiteOnly) for(i in  1:length(datas)) {
      chFini <- is.finite(datas[[i]])
      if(any(!chFini)) datas[[i]] <- datas[[i]][which(chFini)] }
    ## prepare input data
    if(length(col) <n) col <- rep(col,n)[1:n]
    if(horizontal & length(datas) >1) {
      datas <- datas[length(datas):1]                                          # return order that plot will read top -> bottom
      col <- col[length(col):1] }
    if(length(at) != n) at <- 1:n
    su <- sapply(datas,summary)
    if(length(dim(su)) <2) su <- as.matrix(su)
    ra2 <- apply(su[c("Min.","1st Qu.","Median","3rd Qu."),],2, function(x) c(lwr=min(x[1], x[2]-(x[4]-x[2])*rangeVa,na.rm=TRUE), 
      upr=max(x[5], x[4]+(x[4]-x[2])*rangeVa,na.rm=TRUE)) )   # min-max range exteded by 'rangeVa'
    smDensity <- sm::sm.density 
    smArgs <- list(display="none")
    if(!(is.null(hh))) smArgs <- c(smArgs, hh=hh)
    raGlob <- c(min=min(su["Min.", ],na.rm=TRUE), max=max(su["Max.", ],na.rm=TRUE))     # used ?
    smFx <- function(yy,wex,halfViolin) {
      suY <- summary(yy)
    	raY <- c(lwr=min(suY[1],na.rm=TRUE), upr=max(suY[6],na.rm=TRUE) )         # no extension for estimation range
    	smRes <- do.call("smDensity", c(list(yy, xlim=raY), smArgs))
    	## do we care about the upper and lower ends of a variability band, or standard error estimate (which may not always ge produced) ? 
    	smRes$estimate <- smRes$estimate*wex                                       # change proportionally with of violins
      ## truncate to range of real values
    	chRa <- smRes$eval.points > min(suY[1],na.rm=TRUE) & smRes$eval.points < max(suY[6],na.rm=TRUE)
      if(!all(chRa)) {
    	  smRes$eval.points <- smRes$eval.points[which(chRa)]
    	  smRes$estimate <- smRes$estimate[which(chRa)] }
    	## now add points to start & end from baseline=0
    	chBa <- any(smRes$estimate[c(1,length(smRes$estimate))] != 0)
    	if(chBa) { nPo <- length(smRes$estimate)
    	  smRes$estimate <- c(0,smRes$estimate,0)
    	  smRes$eval.points <- smRes$eval.points[c(1,1:nPo,nPo)] }
    	if(!any(sapply(list("yes","pairwise",TRUE),identical,halfViolin ))) { 
        ## need to 'double' data for symmetric violin
        nPo <- length(smRes$estimate)     # otherwise left/upper half
    	  smRes$estimate <- c(smRes$estimate,-1*smRes$estimate[nPo:1]  )
    	  smRes$eval.points <- smRes$eval.points[c(1:nPo,nPo:1)] }
    	list(evalPo=smRes$eval.points,estimate=smRes$estimate,se=smRes$se,summary=suY) }
    vioDat <- lapply(datas,smFx,wex,halfViolin)
    xLim <- c(1 -max(vioDat[[1]]$estimate,na.rm=TRUE), n+max(vioDat[[n]]$estimate,na.rm=TRUE))  # supposed vertical display
    yLim <- if(length(ylim) !=2) range(su[c("Min.","Max."),],na.rm=TRUE) else ylim   	  
    ## configure data/sample-names
    label <- if(is.null(nameSer)) names(datas) else nameSer
    if(is.null(label) | all(is.na(label))) label <- 1:n
    if(identical(halfViolin,"pairwise")) {
      xLim <- xLim -c(0,ceiling(n/2))
      labA <- label[2*(1:(ceiling(n/2))) -1]
      labB <- c(label[2*(1:(floor(n/2)))], if(n %%2 >0) "  ")     
      label <- if(horizontal) paste(labB,labA,sep="\n\n") else paste(labA,labB,sep="    ")
      at <- at[1:ceiling(n/2)] }
    ## configure n per samples
    nDisp <- sapply(datas,function(x) sum(!is.na(x)))
    nDisp <- if(length(unique(nDisp)) >1) paste(rep(c("n=",""), if(n>6) c(1,n-1) else c(n,0)),nDisp,sep="") else paste(" each n=",nDisp[1])
    if(!add) graphics::plot.new()
    ## main plotting 
    if(!horizontal) {
      ## plot vertical
      if(!add) {
        graphics::plot.window(xlim=xLim, ylim=yLim)
        graphics::axis(2,las=las)
        graphics::axis(1, at=at, label=label,las=las)                                                 
        if(!is.null(tit)) graphics::title(main=tit)
      }
      graphics::box()
      if(identical(halfViolin,"pairwise"))  {
        for(i in 1:n) graphics::polygon(ceiling(i/2)+vioDat[[i]]$estimate*(1 -2*(i %% 2)),vioDat[[i]]$evalPo,col=col[i],border=border)
      } else {for(i in 1:n) graphics::polygon(i+vioDat[[i]]$estimate,vioDat[[i]]$evalPo,col=col[i],border=border) }          
    } else {
      ## this is a horizontal plot ..
      if(!add) {
        graphics::plot.window(xlim=yLim, ylim=xLim)
        graphics::axis(1)
        graphics::axis(2, at=at, label=label,las=las)
        if(!is.null(tit)) graphics::title(main=tit)
      }
      graphics::box()
      if(identical(halfViolin,"pairwise"))  {
        for(i in 1:n) graphics::polygon(vioDat[[i]]$evalPo,ceiling(i/2)+vioDat[[i]]$estimate*(1 -2*(i %% 2)),col=col[i],border=border)
      } else {for(i in 1:n) graphics::polygon(vioDat[[i]]$evalPo,i+vioDat[[i]]$estimate,col=col[i],border=border) }          
    }
    ## need to fix 'pairwise' : display of names of data, xLim when horizontal=TRUE
    ## display median : not included yet
    if(length(nDisp) <2) graphics::mtext(nDisp, adj=0, side=3, cex=0.7, line=-1,las=1) else {
      if(identical(halfViolin,"pairwise")) { 
        nDisp <- cbind(nDisp[2*(1:(ceiling(n/2))) -1], c(nDisp[2*(1:(floor(n/2)))],if(n %%2 >0) " "))
        if(horizontal) graphics::mtext(paste(nDisp[,2],nDisp[,1],sep="\n\n"),at=(1:nrow(nDisp))+0.01,side=3-horizontal, cex=0.7,line=-1.5,las=1) else {
          graphics::mtext(paste(nDisp[,1],nDisp[,2],sep=" , "),at=1:nrow(nDisp),side=3, cex=0.7,line=-0.8)}
      } else graphics::mtext(nDisp,at=(1:n)+horizontal/9, side=3-horizontal,las=1, cex=0.7, line=-0.8-horizontal)
    }}
}
    
