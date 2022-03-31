#' Violin-plots version W
#'
#' This function allows generating \href{https://en.wikipedia.org/wiki/Violin_plot}{Violin plots}) using a variety of input formats and offers additional options for colors. 
#' Main input may be multiple vectors, a matrix or list of multiple data-elements (entries may be of variable length),
#'	individual colors for different sets of data or color-gradients can be specified, and the display of n per set of data was integtated 
#' (based on an inspiration from the discussion 'Removing-NAs-from-dataframe-for-use-in-Vioplot' on the forum Nabble). 
#' It is also possible to plot pairwise half-violins for easier pairwise-comparisons (using \code{halfViolin="pairwise"}).
#' Many arguments are kept similar to \href{https://CRAN.R-project.org/package=vioplot}{vioplot} (here, the package \code{vioplot} is not required/used).
#' 
#' @details The (relative) width of the density-profiles ('Violins') may be manually adjusted using the parameter \code{wex} which applieds to all profiles drawn.
#'  Please note that different n (eg for different columns) will not be shown, so far.
#' 
#' Note : Arguments have to be given with full names, lazy evaluation of arguments will not work properly with this function (since '...' is used to capture additional data-sets).    
#' Note : \href{https://CRAN.R-project.org/package=vioplot}{vioplot} offers better options for plotting formulas 
#' 
#' @param x (matrix, list or data.frame) data to plot, or first series of data
#' @param ... (numeric) additional sets of data to plot 
#' @param finiteOnly (logical) eliminate non-finite elements to avoid potential errors (eg when encountering \code{NA})
#' @param halfViolin (logical or character) decide with \code{TRUE} or \code{FALSE} if full or only half of violins should be plotted, if "pairwise" always 2 data-sets will be plotted back-to-back
#' @param boxCol (character) decide if boxplot should be adde inside the violin, use "def" for default transparent grey
#' @param hh (numeric, length <4) smoothing parameter (standard deviation to kernel function, if omited anormal optimal smoothing parameter is used); equivalent to argument \code{h} in package  \href{https://CRAN.R-project.org/package=vioplot}{vioplot} ; see also \code{\link[sm]{sm.density}} 
#' @param ylim (\code{NULL} or numeric, length=2) custom limit on y-axis, see also \code{\link[graphics]{par}} 
#' @param nameSer (character) custom label for data-sets or columns (length must match number of data-sets)
#' @param cexNameSer (numeric) size of individual data-series labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param horizontal (logical) orientation of plot
#' @param col (character or integer) custom colors or gradients like 'rainbow', 'grayscale', 'heat.colors', 'topo.colors', 'Spectral' or 'Paired',  or you may use colors made by the package \href{https://CRAN.R-project.org/package=colorRamps}{colorRamps}  
#' @param border (character) custom color for figure border
#' @param xlab (character) custom x-axis label
#' @param ylab (character) custom y-axis label
#' @param cexLab (numeric) size of axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param cexAxis (numeric) size of numerix axis labels as cex-expansion factor (see also \code{\link[graphics]{par}}) 
#' @param lty (integer) line-type for linear regression line (see also \code{\link[graphics]{par}}) 
#' @param pointCol (character or numeric) display of median: color (defauly white)
#' @param cexPt (numeric) display of median : size of point as cex-expansion factor (see also \code{\link[graphics]{par}}) 
#' @param tit (character) custom title to figure
#' @param las (integer) orientation of axis labels (see also \code{\link[graphics]{par}}) 
#' @param lwd (integer) width of line(s) (see also \code{\link[graphics]{par}}) 
#' @param rectCol (character) color of rectangle
#' @param at (numeric) custom locoation of data-series names, ie the points at which tick-marks are to be drawn, will be passed to \code{\link[graphics]{axis}}, it's length ust match the number of data-sets
#' @param add (logical) add to existing plot if \code{TRUE}
#' @param wex (integer) relative expansion factor of the violin
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @param debug (logical) additional messages for debugging 
#' @return This function plots a figure (to the current graphical device)
#' @seealso the package \href{https://CRAN.R-project.org/package=vioplot}{vioplot}, \code{\link[sm]{sm}} is used for the density estimation
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @examples
#' set.seed(2013)
#' dat6 <- matrix(round(rnorm(300) +3, 1), ncol=6, 
#' 	 dimnames=list(paste0("li",1:50), letters[19:24])) 
#' vioplotW(dat6)
#' ## variable number of elements (each n is displayed)
#' dat6b <- apply(dat6, 2, function(x) x[which(x < 5)])
#' dat6b[[4]] <- dat6b[[4]][dat6b[[4]] < 4]
#' vioplotW(dat6b, col="Spectral")
#' vioplotW(dat6b, col="Spectral" ,halfViolin="pairwise", horizontal=TRUE)
#' vioplotW(dat6b, col="Spectral", halfViolin="pairwise", horizontal=FALSE)
#' @export
vioplotW <- function(x, ..., finiteOnly=TRUE, halfViolin=FALSE, boxCol="def", hh=NULL, ylim=NULL, nameSer=NULL, cexNameSer=NULL, horizontal=FALSE,
	col="rainbow", border="black", xlab=NULL, ylab=NULL, cexLab=NULL, cexAxis=NULL, lty=1, pointCol=NULL, cexPt=NULL,
    tit=NULL, las=1, lwd=1, rectCol="black", at=0, add=FALSE, wex=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## wr variant for violin-plots, inspired&extended based on https://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html
  ## 'col'.. individual fill-colors, may specify gradients (default 'rainbow','grayscale','Spectral','Paired')
  ## default display (short)column names and n
  doPlot <- TRUE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="vioplotW")
  msg <- NULL
  rangeVa <- 1.5                                                   # how much to extent fitted distribution beyond real min/max (for ra2)
  if(is.null(pointCol)) pointCol <- grDevices::rgb(1,1,1,0.9)      # color of (white) point to mark median
  if(is.null(cexPt)) cexPt <- 1.7
  argN <- c(x=deparse(substitute(x)), sup=deparse(substitute(...)))
  chPa <- requireNamespace("sm", quietly=TRUE)
  if(!chPa) { doPlot <- FALSE
    message(fxNa,"Cannot find package 'sm' which is needed for function, please install first from CRAN !")} 
  if(length(x) <1) { warning(fxNa,"Input 'x' seems empty !"); doPlot <- FALSE} 

  if(isTRUE(halfViolin)) halfViolin <- FALSE
  if(!isTRUE(horizontal)) horizontal <- FALSE
  
  fxArg <- c("x","finiteOnly","halfViolin","boxCol","hh","ylim","nameSer","horizontal","col","border","las",
    "lty","tit","lwd","rectCol","at","add","wex","drawRect",  "colMed","pchMed","silent","debug","callFrom") 
  fxAr2 <- union(argN,fxArg[-1])
  chX <- fxAr2 %in% "x"
  if(any(chX)) fxAr2 <- fxAr2[-which(chX)]
  colNx <- if(length(dim(x)) >1) colnames(x) else argN[1]
  if(debug) message(fxNa, "Starting to prepare data, vv1")
  datas <- list(x=x, ...)
  if(debug) {message(fxNa, "Starting to check main input, vv1b"); vv1b <- list(x=x,datas=datas,argN=argN,fxArg=fxArg)}
  datas <- wrMisc::asSepList(datas, asNumeric=TRUE, exclElem=fxAr2, silent=silent, debug=debug, callFrom=fxNa)
  if(length(datas) >0) {
    if(debug) message(fxNa,"vv2, ",length(datas)," additional objects found")
    datas <- lapply(datas, wrMisc::naOmit)
    chLe <- sapply(datas, length)
    if(TRUE & any(chLe <4)) {
      ## simply remove all entries with less than 4 values (makes no sense to plot density
      msg <- paste0(" NOTE: ",sum(chLe <4)," (out of ",length(chLe),") entries have <4 values, eliminated from data to plot")
      if(!silent) message(fxNa, msg)
      datas <- datas[which(chLe >3)]
    }  
  }
  if(length(datas) <1) doPlot <- FALSE
  if(debug) {message(fxNa,"vv3, ready to pepare plotting ? ",doPlot)}
  n <- length(datas)
  if(doPlot & length(n) >0) {  
    if("rainbow" %in% col) col <- grDevices::rainbow(round(n*1.08))[1:n]
    if("grayscale" %in% col) col <- grDevices::gray.colors(n)
    if("heat.colors" %in% col) col <- grDevices::heat.colors(n)
    if("topo.colors" %in% col) col <- grDevices::topo.colors(n)
    
    chPa <- requireNamespace("RColorBrewer", quietly=TRUE)
    if(!chPa & any(c("Spectral","Paired") %in% col)) { col <- "rainbow"
      warning("Cannot find package 'RColorBrewer' (install first from CRAN) which is needed for your choice of argment 'col' !") }
    if("Spectral" %in% col & chPa) col <- RColorBrewer::brewer.pal(min(n,11),"Spectral")
    if("Paired" %in% col & chPa)  col <- RColorBrewer::brewer.pal(12,"Paired")[c(5:6,1:4,7:12)][1:min(n,12)]
    if(length(n) ==1 & length(datas)==1) datas[[1]] <- as.matrix(datas[[1]])
    if(isTRUE(finiteOnly)) for(i in  1:n) {
      chFini <- is.finite(datas[[i]])
      if(any(!chFini)) datas[[i]] <- try(datas[[i]][which(chFini)], silent=TRUE) }
    if(identical(boxCol,"def")) boxCol <- grDevices::rgb(0,0,0,0.4)
    ## prepare input data
    if(length(col) <n) col <- rep(col,n)[1:n]
    if(isTRUE(horizontal) & n >1) {
      datas <- datas[n:1]                                          # return order that plot will read top -> bottom
      col <- col[length(col):1] }
    if(length(at) != n) at <- 1:n
    if(debug) {message(fxNa,"vv4, prepare 'at' ",wrMisc::pasteC(utils::head(at)))}
    su <- if(n==1) as.matrix(summary(as.numeric(datas[[1]]))) else sapply(datas, summary)
    chNum <- is.numeric(su)
    if(!chNum) warning(fxNa," TROUBLE ahead, data from 'x' seem NOT numeric !")
     fxSu <- function(y) c(lwr=min(y[1], y[2] -(y[5] -y[2])*rangeVa, na.rm=TRUE), 
      firQ=as.numeric(y[2]), med=as.numeric(y[3]), thiQ=as.numeric(y[5]), upr=max(y[6], y[5] +(y[5] -y[2])*rangeVa, na.rm=TRUE) ) 
    ra2 <- if(ncol(su) ==1) as.matrix(fxSu(su[,1])) else apply(su, 2, fxSu)
    smDensity <- sm::sm.density 
    smArgs <- list(display="none")
    if(!(is.null(hh))) {smArgs[1 +(1:length(hh))] <- hh; names(smArgs)[length(smArgs) +1 -(length(hh):1)]
      if(debug) message(fxNa, "New arguments added to smArgs :",unlist(hh)) }
    raYGlob <- c(min=min(su["Min.", ], na.rm=TRUE), max=max(su["Max.", ], na.rm=TRUE))     # used ?
    if(debug) message(fxNa,"vv5, raYGlob ",wrMisc::pasteC(signif(raYGlob,4)))

    smFx <- function(yy, halfViolin) {
      if(length(dim(yy)) <2) yy <- as.matrix(yy)
      suY <- if(ncol(yy)==1) summary(as.numeric(yy)) else summary(yy)
    	raY <- c(lwr=min(suY[1], na.rm=TRUE), upr=max(suY[6], na.rm=TRUE) )         # no extension for estimation range
    	smRes <- do.call("smDensity", c(list(yy, xlim=raY), smArgs))
    	## do we care about the upper and lower ends of a variability band, or standard error estimate (which may not always get produced) ? 
      ## truncate to range of real values
    	chRa <- smRes$eval.points >= min(suY[1], na.rm=TRUE) & smRes$eval.points <= max(suY[6], na.rm=TRUE)
      if(!all(chRa)) {
    	  smRes$eval.points <- smRes$eval.points[which(chRa)]
    	  smRes$estimate <- smRes$estimate[which(chRa)] }         
    	## now add points to start & end from baseline=0
    	chBa <- any(smRes$estimate[c(1,length(smRes$estimate))] != 0)
    	if(chBa) { nPo <- length(smRes$estimate)
    	  smRes$estimate <- c(0, smRes$estimate,0)
    	  smRes$eval.points <- smRes$eval.points[c(1,1:nPo,nPo)] }
    	if(!any(sapply(list("yes","pairwise",TRUE), identical, halfViolin ))) { 
        ## need to 'double' data for symmetric violin
        nPo <- length(smRes$estimate)                      # otherwise left/upper half
    	  smRes$estimate <- c(smRes$estimate, -1*smRes$estimate[nPo:1]  )
    	  smRes$eval.points <- smRes$eval.points[c(1:nPo, nPo:1)] }
    	list(evalPo=smRes$eval.points, estimate=smRes$estimate, se=smRes$se, summary=suY) }
     
    vioDat <- lapply(datas, smFx, halfViolin)
    xLim <- c(1 -max(vioDat[[1]]$estimate, na.rm=TRUE), n +max(vioDat[[n]]$estimate, na.rm=TRUE))  # supposed vertical display
    yLim <- if(length(ylim) !=2) range(su[c("Min.","Max."),], na.rm=TRUE) else ylim
    if(debug) message(fxNa,"Init xlim=",wrMisc::pasteC(signif(xLim),4),"  ylim=",wrMisc::pasteC(signif(yLim,4)) )
    
    ## configure data/sample-names
    label <- if(is.null(nameSer)) names(datas) else nameSer
    if(is.null(label) | all(is.na(label))) label <- 1:n
    if(identical(halfViolin,"pairwise")) {
      xLim <- xLim -c(0,ceiling(n/2))
      labA <- label[2*(1:(ceiling(n/2))) -1]
      labB <- c(label[2*(1:(floor(n/2)))], if(n %%2 >0) "  ")     
      ## check if half-series names just differ by extension
      if(FALSE) {
        labelH <- "(half-series names  differing just by extension) not yet implemented"
      } else {      
      ## (if not differing just by extension) concatenate half-series names and  decide if text should contain carriage return      
      nBL <- nBR <- nBl <- nchar(labA) - nchar(labB)
      chLe <- nBL < 0
      if(any(!chLe)) {nBL[which(!chLe)] <- 0       # left text is wider, no need to add add'l blanks)
        nBR[which(!chLe)] <- abs(nBR) }
      if(any(chLe)) {nBL[which(chLe)] <- abs(nBl)  # left text is smaller, need to add extra blanks
        if(any(chLe)) nBR[which(chLe)] <- 0 }
      conMultFx <- function(y) sapply(y, function(x) paste(rep(" ",x),collapse=""))
      label <- if(isTRUE(horizontal)) paste(labB,labA,sep="\n\n") else paste0(conMultFx(nBL),labA, paste0(rep(" ",c(2))), conMultFx(nBR),labB)}   # try to make symmetrix by adding blanks
      at <- at[1:ceiling(n/2)] }
    ## configure n per samples
    nDisp <- sapply(datas, function(x) sum(!is.na(x)))
    nDisp <- if(length(unique(nDisp)) >1) paste0(rep(c("n=",""), if(n >6) c(1,n-1) else c(n,0)),nDisp) else paste(if(n >1) " each","n=",nDisp[1])
    if(!isTRUE(add)) graphics::plot.new()
    ## for boxplot like insert
    boxFa <- 0.06                                  # box half-width
    boxCoor <- cbind(xL=(1:n) -boxFa, yB=ra2[2,], xR=(1:n) +boxFa, yT=ra2[4,])  # basic representation
    ## prepare half-violins
    if(identical(halfViolin,"pairwise")) {
      boxFa <- boxFa +identical(halfViolin,"pairwise")/50              # adjust box width
      pwBoC <- matrix(NA, nrow=5, ncol=n, dimnames=list(c("xLe","yBo","xRi","yTo","xPo"),NULL))  # coords for box
      impa <- which(1:n %% 2 >0)
      pair <- 2* which(impa +1 <= n)
      pwBoC[,impa] <- rbind(1:length(impa) -boxFa, ra2[2,impa], 1:length(impa), ra2[4,impa], 1:length(impa) -boxFa/2)
      pwBoC[,pair] <- rbind(1:length(pair), ra2[2,pair], 1:length(pair) +boxFa, ra2[4,pair], 1:length(impa) +boxFa/2)
    }
    ## automatic adjustement of wex
    ##  make automatic size of wex : smaller if many violins ...
    if(length(wex) !=1) wex <- NULL else if(!is.finite(wex)) wex <- NULL    
    if(length(wex) <1) {     
      we1 <- sapply(vioDat, function(y) stats::quantile(y$estimate, c(0.75,1), na.rm=TRUE))
      wex <- signif(0.5/max(apply(we1, 1, max) *c(0.99)), 4) 
      if(debug) message(fxNa," use wex=",signif(wex,4))
    }
    for(i in 1:length(vioDat)) vioDat[[i]]$estimate <- vioDat[[i]]$estimate *wex      # change proportionally with of violins
    ## main plotting
    if(!isTRUE(horizontal)) {                                # plot vertical
      if(!isTRUE(add)) {
        graphics::plot.window(xlim=xLim, ylim=yLim)
        graphics::axis(2, las=las, cex.axis=cexAxis)        
        graphics::axis(1, at=at, label=label, las=las, cex.axis=cexNameSer, adj=0.5)     # name/labels for indiv series of data
        if(length(xlab) >0) graphics::mtext(xlab, cex=cexLab, side=1, line=2.5) 
        if(length(ylab) >0) graphics::mtext(ylab, cex=cexLab, side=2, line=2.5) 
        if(!is.null(tit)) graphics::title(main=tit)
      }
      graphics::box()
      if(identical(halfViolin,"pairwise"))  {
        for(i in 1:n) graphics::polygon(ceiling(i/2) +vioDat[[i]]$estimate*(1 -2*(i %% 2)), vioDat[[i]]$evalPo, col=col[i], border=border)
        graphics::rect(pwBoC[1,], pwBoC[2,], pwBoC[3,], pwBoC[4,], col=boxCol, border=NA)
        graphics::points(pwBoC[5,],ra2[3,], pch=16, cex=cexPt, col=pointCol)
        graphics::segments(pair/2, ra2[2,impa], pair/2, ra2[4,impa], col="grey", lty=2)     # grey line to separate adjacent boxes
      } else {for(i in 1:n) graphics::polygon(i+vioDat[[i]]$estimate, vioDat[[i]]$evalPo, col=col[i], border=border) 
        graphics::rect(boxCoor[,1], boxCoor[,2], boxCoor[,3], boxCoor[,4], col=boxCol, border=NA)
        graphics::points(1:n, ra2[3,], pch=16, cex=cexPt, col=pointCol) 
      }          
      if(length(boxCol >0)) {                      ## add boxplot-like
      }
    } else {
      ## this is a horizontal plot ..
      if(!isTRUE(add)) {
        graphics::plot.window(xlim=yLim, ylim=xLim)
        graphics::axis(1, cex.axis=cexNameSer)
        graphics::axis(2, cex.axis=cexAxis, at=at, label=label,  cex=cexNameSer, las=las)
        if(length(xlab) >0) graphics::mtext(xlab, cex=cexLab, side=1, line=2.5) 
        if(length(ylab) >0) graphics::mtext(ylab, cex=cexLab, side=2, line=2.5) 
        if(!is.null(tit)) graphics::title(main=tit)
      }
      graphics::box()
      if(identical(halfViolin,"pairwise"))  {
        for(i in 1:n) graphics::polygon(vioDat[[i]]$evalPo, ceiling(i/2) +vioDat[[i]]$estimate*(1 -2*(i %% 2)), col=col[i], border=border)
        graphics::rect(pwBoC[2,], pwBoC[1,], pwBoC[4,], pwBoC[3,], col=boxCol, border=NULL)
        graphics::points(ra2[3,], pwBoC[5,], pch=16, cex=cexPt, col=pointCol) 
        graphics::segments(ra2[2,impa], pair/2, ra2[4,impa], pair/2, col="grey",lty=2)     # grey line to separate adjacent boxes
      } else {for(i in 1:n) graphics::polygon(vioDat[[i]]$evalPo,i +vioDat[[i]]$estimate, col=col[i], border=border) 
        graphics::rect(boxCoor[,2], boxCoor[,1], boxCoor[,4], boxCoor[,3], col=boxCol, border=NA)
        graphics::points( ra2[3,], 1:n, pch=16, cex=cexPt, col=pointCol) 
        }
          
    }
    ## display n by series (only if variable along data-set)
    if(length(nDisp) <2) graphics::mtext(nDisp, adj=0, side=3, cex=0.7, line=-1, las=1) else {
      if(identical(halfViolin,"pairwise")) { 
        nDisp <- cbind(nDisp[2*(1:(ceiling(n/2))) -1], c(nDisp[2*(1:(floor(n/2)))], if(n %%2 >0) " "))
        if(horizontal) graphics::mtext(paste(nDisp[,2],nDisp[,1],sep="\n\n"), at=(1:nrow(nDisp))+0.01, side=3-horizontal, cex=0.7, line=-1.5, las=1) else {
          graphics::mtext(paste(nDisp[,1],nDisp[,2],sep=" , "), at=1:nrow(nDisp), side=3, cex=0.7, line=-0.8)}
      } else graphics::mtext(nDisp, at=(1:n)+horizontal/9, side=3-horizontal, las=1, cex=0.7, line=-0.8-horizontal)
    }
  if(length(msg) >0) graphics::mtext(msg, side=1, cex=0.7,las=1)    # add note about groups/columns of data not plotted
  }
}
    
