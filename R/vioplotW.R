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
#' @param removeEmpty (logical) omit empty series (or less than 4 finite numeric entries) of data from plot
#' @param halfViolin (logical or character) decide with \code{TRUE} or \code{FALSE} if full or only half of violins should be plotted, if "pairwise" always 2 data-sets will be plotted back-to-back
#' @param boxCol (character) decide if boxplot should be adde inside the violin, use "def" for default transparent grey
#' @param hh (numeric, length <4) smoothing parameter (standard deviation to kernel function, if omited anormal optimal smoothing parameter is used); equivalent to argument \code{h} in package  \href{https://CRAN.R-project.org/package=vioplot}{vioplot} ; see also \code{\link[sm]{sm.density}}
#' @param xlim (\code{NULL} or numeric, length=2) custom limit on x-axis, see also \code{\link[graphics]{par}}
#' @param ylim (\code{NULL} or numeric, length=2) custom limit on y-axis, see also \code{\link[graphics]{par}}
#' @param nameSer (character) custom label for data-sets or columns (length must match number of data-sets)
#' @param cexNameSer (numeric) size of individual data-series labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param horizontal (logical) orientation of plot
#' @param col (character or integer) custom colors or gradients like 'rainbow', 'grayscale', 'heat.colors', 'topo.colors', 'Spectral' or 'Paired',  or you may use colors made by the package \href{https://CRAN.R-project.org/package=colorRamps}{colorRamps}
#' @param border (character) custom color for figure border
#' @param xlab (character) custom x-axis label
#' @param ylab (character) custom y-axis label
#' @param cexLab (numeric) size of axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param cexAxis (numeric) size of numeric y-axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})
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
vioplotW <- function(x, ..., finiteOnly=TRUE, removeEmpty=FALSE, halfViolin=FALSE, boxCol="def", hh=NULL, xlim=NULL, ylim=NULL, nameSer=NULL, cexNameSer=NULL, horizontal=FALSE,
	col="rainbow", border="black", xlab=NULL, ylab=NULL, cexLab=NULL, cexAxis=NULL, lty=1, pointCol=NULL, cexPt=NULL,
    tit=NULL, las=1, lwd=1, rectCol="black", at=0, add=FALSE, wex=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## wr variant for violin-plots, inspired&extended based on https://r.789695.n4.nabble.com/Removing-NAs-from-dataframe-for-use-in-Vioplot-td4720274.html
  ## 'col'.. individual fill-colors, may specify gradients (default 'rainbow','grayscale','Spectral','Paired')
  ## default display (short)column names and n
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="vioplotW")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE

  fixArg <- c("x","finiteOnly","halfViolin","boxCol","hh","xlim","ylim","nameSer","horizontal","col","border","las",          # fixed argument names to check (and adjust)
    "lty","tit","lwd","rectCol","at","add","wex","drawRect",  "colMed","pchMed","silent","debug","callFrom")
  datas <- list(...)
       out <- useColNo <- NULL
  argN <- c(x=deparse(substitute(x)), sup=deparse(substitute(...)))
  argL <- match.call(expand.dots = FALSE)$...            # extr arg names, based on https://stackoverflow.com/questions/55019441/deparse-substitute-with-three-dots-arguments
  if(debug) {message(callFrom," -> vioplotW  vv0"); vv0 <- list(datas=datas,x=x,argN=argN,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argL=argL,fixArg=fixArg)}

  ## separate specific arguments from all-input (lazy fitting) & clean datas
  if(length(datas) >0) {
    pMa <- pmatch(names(argL), fixArg)                     # will only find if unique (partial) match
    if(length(argL) >0 && length(pMa) >0) { pMa[which(nchar(names(argL)) <3)] <- NA                # limit to min 3 chars length for lazy matching
      if(any(!is.na(pMa))) for(i in which(!is.na(pMa))) {assign(fixArg[pMa[i]], datas[[i]]); names(datas)[i] <- "replaceReplace"}}
    chRepl <- names(datas) %in% "replaceReplace"
    if(any(chRepl)) {datas <- datas[-which(chRepl)]; argL <- argL[-which(chRepl)]}
  }

  ## combine datas and x & set x to beginning , convert matrices to list of vectors
  if(length(x) >0)  {
    ## check x
    x <- try(wrMisc::asSepList(x, asNumeric=TRUE, silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
    if(inherits(x, "try-error")) {datas <- NULL; if(!silent) message(fxNa,"Problem with separating data from main entry '",argN[1],"'")}}

  if(length(datas) >0)  {
    datas <- try(wrMisc::asSepList(datas, asNumeric=TRUE, silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
    if(inherits(datas, "try-error")) {datas <- NULL; if(!silent) message(fxNa,"Problem with separating data from supplemental data")}}

  if(length(datas) >0)  {
    ## now add x
    if(length(x) >0) {
      datas[[length(datas) +(1:length(x))]] <- x
      newNa <- if(length(argN[1]) >0) argN[1] else "x"
      names(datas)[1 +length(datas) -(length(x):1)] <- if(length(x)==1) newNa else { paste0(newNa,"_",1:length(x))}
      if(length(datas) >1) datas <- datas[c(1 +length(datas) -length(x):1, 2:(length(datas) -length(x)))] }   # reset order
  } else { datas <- x
    if(debug) {message(fxNa,"vv1b")} }

  if(debug) {message(fxNa,"vv1c"); vv1c <- list(datas=datas,x=x,argN=argN,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,removeEmpty=removeEmpty,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argL=argL,fixArg=fixArg)}

  .filtX <- function(x, nMin=4, nullRet=NULL) {x <- wrMisc::naOmit(x); if(length(x) >= nMin) {chN <- is.finite(x); if(sum(chN) >= nMin) x[which(chN)] else nullRet } else nullRet}
  doPlot <- length(x) >0 && length(datas) >0
  if(doPlot) {    ## add'l filter
    ## filter for min length=4
    iniLe <- length(datas)
    datas <- lapply(datas, .filtX, nMin=4, nullRet=NULL)
    datas <- lapply(datas, function(x) {x <- wrMisc::naOmit(x); if(length(x) >3) x})
    if(isTRUE(removeEmpty)) {
      chL <- sapply(datas, length) <1
      if(any(chL)) { if(all(chL)) { warning(fxNa,"INSUFFICIENT data, can't plot"); doPlot <- FALSE
        } else {datas <- datas[which(!chL)]; if(!silent) message(fxNa,"REMOVING ",sum(!chL)," out of ",iniLe, " entries since too short for plotting !")}}
    }
  }


  ##
  msg <- NULL
  rangeVa <- 1.5                                                   # how much to extent fitted distribution beyond real min/max (for ra2)
  if(is.null(pointCol)) pointCol <- grDevices::rgb(1,1,1,0.9)      # color of (white) point to mark median
  if(is.null(cexPt)) cexPt <- 1.7
  pwBoC <- boxCoor <- NULL

  chPa <- requireNamespace("sm", quietly=TRUE)
  if(!chPa) { doPlot <- FALSE
    message(fxNa,"Cannot find package 'sm' which is needed for function, please install first from CRAN !")}
  if(isTRUE(halfViolin)) halfViolin <- FALSE
  if(!isTRUE(horizontal)) horizontal <- FALSE
  if(debug) {message(fxNa,"Starting to prepare data, vv2"); vv2 <- list(datas=datas,x=x,argN=argN,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argN=argN,argL=argL,halfViolin=halfViolin,horizontal=horizontal)}

  n <- length(datas)
  if(doPlot && n >0) {
    if("rainbow" %in% col) col <- grDevices::rainbow(round(n*1.08))[1:n]
    if("grayscale" %in% col) col <- grDevices::gray.colors(n)
    if("heat.colors" %in% col) col <- grDevices::heat.colors(n)
    if("topo.colors" %in% col) col <- grDevices::topo.colors(n)

    chPa <- requireNamespace("RColorBrewer", quietly=TRUE)
    if(!chPa && any(c("Spectral","Paired") %in% col)) { col <- "rainbow"
      warning("Cannot find package 'RColorBrewer' (install first from CRAN) which is needed for your choice of argument 'col' ! (ignoring") }
    if("Spectral" %in% col && chPa) col <- RColorBrewer::brewer.pal(min(n,11), "Spectral")
    if("Paired" %in% col && chPa)  col <- RColorBrewer::brewer.pal(12,"Paired")[c(5:6,1:4,7:12)][1:min(n,12)]

    ## need to re-write ??
    if(length(n) ==1 && n==1) datas[[1]] <- as.matrix(datas[[1]])
    if(isTRUE(finiteOnly)) datas <- lapply(datas, function(x) {chFini <- is.finite(x); if(any(!chFini)) x <- x[which(chFini)]; x})
    if(debug) {message(fxNa, "Starting to prepare data, vv3")}

    if(identical(boxCol,"def")) boxCol <- grDevices::rgb(0,0,0,0.4)
    ## prepare input data
    if(length(col) <n) col <- rep(col,n)[1:n]
    if(isTRUE(horizontal) && n >1) {
      datas <- datas[n:1]                                          # return order that plot will read top -> bottom
      col <- col[length(col):1] }
    if(length(at) != n) at <- 1:n
    if(length(datas)==1) at <- 1
    if(debug) {message(fxNa,"vv5, prepare 'at' ",wrMisc::pasteC(utils::head(at))); vv5 <- list(datas=datas,x=x,argN=argN,rangeVa=rangeVa,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,col=col,at=at,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argN=argN,argL=argL,halfViolin=halfViolin,horizontal=horizontal,n=n)}
    sumry <- function(x) {     # take summary and then points for boxplot  lwr,firQ,med,thiQ,upr
      outNA <- c(lwr=NA,firQ=NA,med=NA,thiQ=NA,upr=NA)
      if(length(x) >0) { if(is.numeric(x)) {y <- summary(as.numeric(x)); c(lwr=min(y[1], y[2] -(y[5] -y[2])*rangeVa, na.rm=TRUE),
      firQ=as.numeric(y[2]), med=as.numeric(y[3]), thiQ=as.numeric(y[5]), upr=max(y[6], y[5] +(y[5] -y[2])*rangeVa, na.rm=TRUE))} else outNA} else outNA
    }
    ra2 <- as.matrix(sapply(datas, sumry))
    colnames(ra2) <- names(datas)

    if(debug) {message(fxNa,"vv5b, prepare 'at' ",wrMisc::pasteC(utils::head(at))); vv5b <- list()}
    smDensity <- sm::sm.density
    smArgs <- list(display="none")
    if(!(is.null(hh))) {smArgs[1 +(1:length(hh))] <- hh; names(smArgs)[length(smArgs) +1 -(length(hh):1)]
      if(debug) message(fxNa, "New arguments added to smArgs :",unlist(hh)) }
    raYGlob <- c(min=min(ra2[1,], na.rm=TRUE), max=max(ra2[5,], na.rm=TRUE))     # used ?
    if(debug) {message(fxNa,"vv6, raYGlob ",wrMisc::pasteC(signif(raYGlob,4)))
      vv6 <- list(datas=datas,x=x,argN=argN,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,smDensity=smDensity,smArgs=smArgs,raYGlob=raYGlob,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argN=argN,argL=argL,halfViolin=halfViolin,horizontal=horizontal)}

    smFx <- function(yy, halfViolin, debug=FALSE) {
      if(length(dim(yy)) <2) yy <- as.matrix(yy)
      suY <- if(ncol(yy)==1) summary(as.numeric(yy)) else summary(yy)
    	raY <- c(lwr=min(suY[1], na.rm=TRUE), upr=max(suY[6], na.rm=TRUE) )         # no extension for estimation range
    	smRes <- do.call("smDensity", c(list(yy, xlim=raY), smArgs))
        if(debug) {message("smFx : sm1"); sm1 <- list(yy=yy,halfViolin=halfViolin,suY=suY,raY=raY,smRes=smRes)}
    	## do we care about the upper and lower ends of a variability band, or standard error estimate (which may not always get produced) ?
      ## truncate to range of real values
    	chRa <- smRes$eval.points >= min(suY[1], na.rm=TRUE) & smRes$eval.points <= max(suY[6], na.rm=TRUE)
      if(!all(chRa)) {
    	  smRes$eval.points <- smRes$eval.points[which(chRa)]
    	  smRes$estimate <- smRes$estimate[which(chRa)] }
        if(debug) {message("smFx :sm2"); sm2 <- list(yy=yy,halfViolin=halfViolin,suY=suY,raY=raY,smRes=smRes)}
    	## now add points to start & end from baseline=0
    	chBa <- any(smRes$estimate[c(1,length(smRes$estimate))] != 0)
    	if(chBa) { nPo <- length(smRes$estimate)
    	  smRes$estimate <- c(0, smRes$estimate,0)
    	  smRes$eval.points <- smRes$eval.points[c(1,1:nPo,nPo)] }
        if(debug) {message("smFx : sm3"); sm3 <- list(yy=yy,halfViolin=halfViolin,suY=suY,raY=raY,smRes=smRes)}
    	if(!any(sapply(list("yes","pairwise",TRUE), identical, halfViolin ))) {
        ## need to 'double' data for symmetric violin
        nPo <- length(smRes$estimate)                      # otherwise left/upper half
    	  smRes$estimate <- c(smRes$estimate, -1*smRes$estimate[nPo:1]  )
    	  smRes$eval.points <- smRes$eval.points[c(1:nPo, nPo:1)] }
    	list(evalPo=smRes$eval.points, estimate=smRes$estimate, se=smRes$se, summary=suY) }

    smFx2 <- function(yy, halfViolin, debug=FALSE, callFrom=NULL) {
      out <- try(smFx(yy, halfViolin, debug=FALSE), silent=TRUE)
      #chEr <- sapply(out, inherits, "try-error")
      #if(any(chEr)) { warning(callFrom,"Problem with calculating density ! Need to omit this series")
      #  out <- out[which(!chEr)] }
      out }

    vioDat <- lapply(datas, smFx2, halfViolin, debug=debug, callFrom=fxNa)
    ## correct for empty series
    chNull <- colSums(is.na(ra2)) ==5
    if(any(chNull)) for(i in which(chNull)) vioDat[[i]] <- list(estimate=NULL)

    xLim <- if(length(xlim) ==2) xlim else {               # 
      c(1 -max(0.2, vioDat[[1]]$estimate, na.rm=TRUE), n +max(0.2, vioDat[[n]]$estimate, na.rm=TRUE))}  # supposed vertical display, note final xLim is changed by wex
    yLim <- if(length(ylim) ==2) ylim else range(ra2[c(1,5),], na.rm=TRUE)
    if(debug) {message(fxNa,"vv7, Init xlim=",wrMisc::pasteC(signif(xLim),4),"  ylim=",wrMisc::pasteC(signif(yLim,4)) )
      vv7 <- list(datas=datas,x=x,vioDat=vioDat,xLim=xLim,yLim=yLim,argN=argN,finiteOnly=finiteOnly,halfViolin=halfViolin,boxCol=boxCol,hh=hh,smDensity=smDensity,smArgs=smArgs,raYGlob=raYGlob,xlim=xlim,ylim=ylim,nameSer=nameSer,tit=tit,argN=argN,argL=argL,halfViolin=halfViolin,horizontal=horizontal)}

    ## configure data/sample-names
    label <- if(is.null(nameSer)) names(datas) else nameSer
    if(is.null(label) || all(is.na(label))) label <- 1:n
    if(identical(halfViolin,"pairwise")) {
      xLim <- xLim -c(0, ceiling(n/2))
      labA <- label[2*(1:(ceiling(n/2))) -1]
      labB <- c(label[2*(1:(floor(n/2)))], if(n %%2 >0) "  ")
      ## check if half-series names just differ by extension
      if(FALSE) {
        labelH <- "(half-series names  differing just by extension) - not yet implemented"
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
     if(debug) {message(fxNa,"vv7c"); vv7c <- list( )}

    ## automatic adjustement of wex
    ##  make automatic size of wex : smaller if many violins ...
    if(length(wex) !=1) wex <- NULL else if(!is.finite(wex)) wex <- NULL
    if(length(wex) <1) {
      we1 <- sapply(vioDat, function(y) if(length(y$estimate) >0) stats::quantile(y$estimate, c(0.75,1), na.rm=TRUE) else rep(NA,2))
      wex <- if(length(dim(we1)) >1) signif(0.5/max(apply(we1, 1, max, na.rm=TRUE) *c(0.99), na.rm=TRUE), 4) else signif(0.5/max(we1, na.rm=TRUE))
      if(debug) {message(fxNa,"vv7d  use wex=",signif(wex,4))}
    }
    for(i in 1:length(vioDat)) vioDat[[i]]$estimate <- vioDat[[i]]$estimate *wex      # change proportionally with of violins
    if(debug) {message(fxNa,"vv8"); vv8 <- list(datas=datas, vioDat=vioDat,ra2=ra2,pwBoC=pwBoC,boxCoor=boxCoor,xLim=xLim,ylim=yLim,at=at,wex=wex,halfViolin=halfViolin,cexPt=cexPt,boxCol=boxCol )}
    ## main plotting
    if(!isTRUE(horizontal)) {                                # plot vertical
      xLim <- c(1- max(vioDat[[i]]$estimate)*1.1, n+ max(vioDat[[n]]$estimate)*1.1)              # readjust xLim to final wex
      cat("..new xLim ",round(xLim,3),"\n")
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
      } else {for(i in 1:n) graphics::polygon(i +vioDat[[i]]$estimate, vioDat[[i]]$evalPo, col=col[i], border=border)
        graphics::rect(boxCoor[,1], boxCoor[,2], boxCoor[,3], boxCoor[,4], col=boxCol, border=NA)
        graphics::points(1:n, ra2[3,], pch=16, cex=cexPt, col=pointCol)
      }
      if(length(boxCol >0)) {                      ## add boxplot-like
      }
    } else {
      ## this is a horizontal plot ..
      yLim <- c(1- max(vioDat[[i]]$estimate)*1.1, n+ max(vioDat[[n]]$estimate)*1.1)               # readjust xLim to final wex
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
      } else graphics::mtext(nDisp, at=(1:n) +horizontal/9, side=3-horizontal, las=1, cex=0.7, line=-0.8-horizontal)
    }
  if(length(msg) >0) graphics::mtext(msg, side=1, cex=0.7, las=1)    # add note about groups/columns of data not plotted
  }
}
    
