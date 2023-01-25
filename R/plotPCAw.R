#' PCA plot with bag-plot to highlight groups
#'
#' @description
#' This function allows to plot \href{https://en.wikipedia.org/wiki/Principal_component_analysis}{principal components analysis (PCA)},
#' with options to show center and potential outliers for each of the groups (columns of data).
#' The main points of this implementation consist in offering bagplots to highlight groups of columns/samples and support to (object-oriented)
#'  output from \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} and \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo}.
#'
#' @details
#' One motivation for this implementation of plotting PCA was to provide a convenient way for doing so with of MArrayLM-objects or lists
#' as created by \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} and \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo}.
#'
#' Another motivation for this implementation come from integrating the idea of bag-plots to better visualize different groups of points
#' (if they can be organized so beforehand as distinct groups) :
#' The main body of data is shown as 'bag-plots' (a bivariate boxplot, see \href{https://en.wikipedia.org/wiki/Bagplot}{Bagplot})
#' with different transparent colors to highlight the core part of different groups (if they contain more than 2 values per group).
#' Furthermore, group centers are shown as average or median (see 'nGrpForMedian') with stars & index-number (if <25 groups).
#'
#' Layout is automatically set to 2 or 4 subplots (if plotting more than 2 principal components makes sense).
#'
#' Note : This function uses \code{\link[stats]{prcomp}} for calculating Eigenvectors and principal components, with default \code{center=TRUE} and \code{scale.=FALSE} (different to \code{princomp()}. which standardizes by default).
#' This way the user has to option to intervene on arguments \code{center} and \code{scale.}. However, this should be done with care.
#'
#' Note: \code{NA}-values cannot (by definition) be processed by (any) PCA - all lines with any non-finite values/content (eg \code{NA}) will be omitted !
#'
#' Note : Package RColorBrewer may be used if available.
#'
#' For more options with PCA (and related methods) you may also see also the package  \href{https://CRAN.R-project.org/package=FactoMineR}{FactoMineR}
#' which provides a very wide spectrum of possibiities, in particular for combined numeric and categorical data.
#'
#' @param dat (matrix, data.frame, MArrayLM-object or list) data to plot. Note: \code{NA}-values cannot be processed - all lines with non-finite data (eg \code{NA}) will be omitted !
#'  In case of  MArrayLM-object or list \code{dat} must conatain list-element named 'datImp','dat' or 'data'.
#' @param sampleGrp (character or factor) should be factor describing groups of replicates, NAs are not supported
#' @param tit (character) custom title
#' @param useSymb (integer) symbols to use (see also \code{\link[graphics]{par}})
#' @param center (logical or numeric) decide if variables should be shifted to be zero centered, argument passed to \code{\link[stats]{prcomp}}
#' @param scale. (logical or numeric) decide if scaling to obtain unit variance, argument passed to \code{\link[stats]{prcomp}}
#'  Alternatively, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.
#' @param colBase (character or integer) use custom colors
#' @param useSymb2 (integer) symbol to mark group-center (no mark of group-center if default NULL) (equivalent to \code{pch}, see also \code{\link[graphics]{par}})
#' @param cexTxt (integer) expansion factor for text (see also \code{\link[graphics]{par}})
#' @param cexSub (integer) expansion factor for subtitle line text (see also \code{\link[graphics]{par}})
#' @param displBagPl (logical) if \code{TRUE}, show bagPlot (group-center) if >3 points per group otherwise the average-confidence-interval
#' @param outCoef (numeric) parameter for defining outliers, see  \code{\link{addBagPlot}} (equivalent to \code{range} in \code{\link[graphics]{boxplot}})
#' @param getOutL (logical) return outlyer samples/values
#' @param showLegend (logical or character) toggle to display legend, if character it designes the location within the plot to display the legend ('bottomleft','topright', etc..)
#' @param nGrpForMedian (integer) decide if group center should be displayed via its average or median value: If group has less than 'nGrpForMedian' values, the average will be used, otherwise the median; if \code{NULL} no group centers will be displayed
#' @param pointLabelPar (character) define formatting for optional labels next to points in main figure (ie PC1 vs PC2); may be \code{TRUE} or list containing elments 'textLabel', 'textCol', 'textCex',
#'  'textOffSet', 'textAdj' for fine-tuning
#' @param rowTyName (character) for subtitle : specify nature of rows (genes, proteins, probesets,...)
#' @param rotatePC (integer) optional rotation (by -1) for fig&ure of the principal components specified by index
#' @param suplFig (logical) to include plots vs 3rd principal component (PC) and Screeplot
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function make a plot and may retiurn an optional matrix of outlyer-data (depending on argument \code{getOutL})
#' @seealso \code{\link[stats]{prcomp}} (used here for the PCA underneith) , \code{\link[stats]{princomp}}, see the package \href{https://CRAN.R-project.org/package=FactoMineR}{FactoMineR} for multiple plotting options or ways of combining categorical and numeric data
#' @examples
#' set.seed(2019); dat1 <- matrix(round(c(rnorm(1000), runif(1000,-0.9,0.9)),2),
#'   ncol=20, byrow=TRUE) + matrix(rep(rep(1:5,6:2), each=100), ncol=20)
#' biplot(prcomp(dat1))        # traditional plot
#' (grp = factor(rep(LETTERS[5:1],6:2)))
#' plotPCAw(dat1, grp)
#' @export
plotPCAw <- function(dat, sampleGrp, tit=NULL, useSymb=c(21:25,9:12,3:4), center=TRUE, scale.=TRUE, colBase=NULL, useSymb2=NULL,
  displBagPl=TRUE, outCoef=2, getOutL=FALSE, cexTxt=1, cexSub=0.6, showLegend=TRUE, nGrpForMedian=6, pointLabelPar=NULL, rowTyName="genes",
  rotatePC=NULL, suplFig=TRUE, callFrom=NULL, silent=FALSE, debug=FALSE) {
  ## note : so far not well adopted for plotting all points in transparent grey with same symb for interactive mouse-over
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotPCAw")
  msg <- " 'dat' must be matrix, data.frame (with at least 1 line and 2 columns) or list/MArrayLM-object (with element 'datImp','dat' or 'data' as matrix or data.frame)"
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if("MArrayLM" %in% class(dat) | is.list(dat)) { chDat <- c("datImp","dat","data","quant","datRaw")  %in% names(dat)
    chGr <- c("groups","grp","sampleGroups")  %in% names(dat)
    if(any(chGr, na.rm=TRUE) & length(sampleGrp) <1) {
      chGrpNa <- c("groups","grp","sampleGroups")
      sampleGrp <- dat[[which(names(dat) %in% chGrpNa)[1]]]
    } else if("sampleSetup" %in% names(dat)) if(any(chGrpNa %in% names(dat$sampleSetup), na.rm=TRUE)) {
      sampleGrp <- dat$sampleSetup[[ which(chGrpNa %in% names(dat$sampleSetup))[1]]] }
    if(any(chDat, na.rm=TRUE)) dat <- dat[[which(names(dat) %in% c("datImp","dat","data","datRaw"))[1]]] else stop(msg) }
  if(any(length(dim(dat)) !=2, dim(dat) < 1:2, na.rm=TRUE)) stop(msg)
  if(length(sampleGrp) <1 | sum(is.na(sampleGrp))==length(sampleGrp)) stop("Argument 'sampleGrp' seems to be empty (or all NA)")
  if(length(sampleGrp) != ncol(dat)) stop(fxNa,"Length of argument 'sampleGrp' must match number of columns of 'dat'")
  if(!is.factor(sampleGrp)) sampleGrp <- as.factor(sampleGrp)
  if(!isFALSE(center)) center <- TRUE
  if(!isFALSE(scale.)) scale. <- TRUE
  if(!isFALSE(suplFig)) suplFig <- TRUE

  nGrpLim <- 25                      # decide/limit if group-centers are shown as stars
  ## reduce to columns defined by 'sampleGrp' (groups declared as NA will be removed) and to columns with finite data
  if(length(dat) >0) {
    dat <- dat[,!is.na(sampleGrp)]
    chFin <- is.finite(dat)
    chCol <- colSums(chFin) <2
    if(sum(chCol) >0){
      dimIni <- dim(dat)
      if(!silent) message(fxNa,"Removing ",sum(chCol)," (out of ",ncol(dat),") columns with less than 2 finite values")
      dat <- dat[,which(!chCol)]
      chFin <- chFin[,which(!chCol)]
      if(ncol(dat) <2) stop("After removing bad (non-finite) columns less than 2 columns remain !")
      sampleGrp <- sampleGrp[which(!chCol)]
      if(length(useSymb) == dimIni[2]) useSymb <- useSymb[which(!chCol)] }
    opar <- graphics::par(no.readonly=TRUE)
    if(suplFig & ncol(dat) > 3) on.exit(opar)  else {
      opar2 <- opar
      opar2 <- opar[which(names(opar) %in% c("mar","mfrow"))]
      on.exit(opar2)
    }
    if(debug) {message(fxNa,"plotPCAw_1") }
    ## remove rows with NAs
    chRow <- rowSums(!chFin)                      # check for lines with non-finite values (avoid error during svd/pca ..)
    if(any(chRow >0, na.rm=TRUE)) if(sum(chRow >0) ==nrow(dat)) stop("All lines have some non-finite values !?!") else {
      if(!silent) message(fxNa,"Eliminating ",sum(chRow >0)," (out of ",nrow(dat),") lines with non-finite values")
      dat <- dat[which(chRow <1),] }
    ## remove 0 variance rows
    isZeroVar <- function(x) length(unique(wrMisc::naOmit(x)))==1
    ch1 <- apply(dat, 1, isZeroVar)
    if(any(ch1, na.rm=TRUE)) {
      if(!silent) message(fxNa,"Eliminating ",sum(ch1)," (out of ",nrow(dat),") lines with 0 variance values (ie all values identical)")
      dat <- dat[which(ch1),] }

    averNa3 <- wrMisc::naOmit(sampleGrp)
    if(suplFig & ncol(dat) > 3) graphics::layout(matrix(c(1,1:3,4,4),3,2,byrow=TRUE), heights=c(10,5,4))
    if(!identical(ncol(dat),length(wrMisc::naOmit(sampleGrp))) & !silent) warning("'plotPCA': number of columns and number of (QC-filtered) arrays in 'AverNa2' don't match !")
    nGrp <- length(levels(sampleGrp))
    if(!is.null(colBase)) if(length(colBase) != nGrp) {
      if(!silent) message(fxNa,"Ignore invalid 'colBase', expecting length ",nGrp)
      colBase <- NULL }
    chPa <- requireNamespace("RColorBrewer", quietly=TRUE)
    if(!silent & !chPa) message(fxNa,"Please install package 'RColorBrewer' for improved color-gradients")
    if(is.null(colBase)) colBase <- if(nGrp <10 & chPa) { RColorBrewer::brewer.pal(max(3,nGrp),"Set1")[1:nGrp]
      } else {
      grDevices::rainbow(1.5 +nGrp*1.08)[1:(nGrp+1)][-ceiling(length(levels(sampleGrp))/2)]}
    useCol <- colBase[order(as.numeric(unique(averNa3)))]
    chGr <- duplicated(sampleGrp)
    if(!any(chGr, na.rm=TRUE)) { displBagPl <- FALSE; if(debug) message(fxNa,"All samples belong to different groups, can't plot bagplots ...")}
    if(!identical(unique(table(averNa3)), as.integer(1))) useCol <- useCol[as.numeric(wrMisc::naOmit(sampleGrp))]
    if(debug) {message(fxNa,"plotPCAw_3") }
  }


  ## MAIN
  if(length(dat >0)) { pca <- try(stats::prcomp(t(dat), center=center, scale.=scale.), silent=TRUE)
    if(inherits(pca, "try-error")) {dat <- NULL; pca <- NULL}}

  if(length(dat) >0) {
    percVar <- 100*round(pca$sdev^2/sum(pca$sdev^2),3)
    pca.grpMed <- pca.grpAv <- matrix(nrow=nGrp, ncol=min(4,ncol(pca$x)))
    tabAverNa3 <- table(averNa3)
    if(length(nGrpForMedian) >0) { for(i in 1:min(4,ncol(pca$x))) {
        pca.grpMed[,i] <- unlist(tapply(pca$x[,i], averNa3, stats::median, na.rm=TRUE))
        pca.grpAv[,i] <- unlist(tapply(pca$x[,i], averNa3, mean, na.rm=TRUE)) }
      rownames(pca.grpMed) <- rownames(pca.grpAv) <- unique(averNa3)
      pca.grpMed[names(tabAverNa3)[which(tabAverNa3 <nGrpForMedian)],] <- pca.grpAv[names(tabAverNa3)[which(tabAverNa3 < nGrpForMedian)],]}
    useCex <- round( (1/nGrp +0.65)*cexTxt, 2)
    useSymb.ori <- rep(useSymb, 1 +length(levels(sampleGrp)) %/% length(useSymb)) [1:length(levels(sampleGrp))] #
    useSymb <- (useSymb.ori[order(unique(as.numeric(averNa3)))] )[as.numeric(wrMisc::naOmit(sampleGrp))]
    useTxtCol <- 1
    figNumOffs <- function(pca, cols=c(1,2), fact=100) signif(c(abs(diff(range(pca$x[,cols[1]])))/fact, abs(diff(range(pca$x[,cols[2]])))/fact), digits=3)
    outL <- list()
    length(outL) <- length(unique(averNa3)); names(outL) <- unique(averNa3)
    lab123 <- paste0("PC",1:3," (",percVar[1:3],"%)")
    useTit <- if(is.null(tit)) "Principal Components of Samples" else tit
    if(suplFig) useTit <- paste(useTit, if(nchar(useTit) < 35) ": 1st and 2nd Component" else ", PC1 & PC2")
    ## optional rotate axis
    if(length(rotatePC) >0 & is.numeric(rotatePC)) {
      chRo <- rotatePC >0 & rotatePC <= ncol(dat)
      if(any(!rotatePC, na.rm=TRUE)) rotatePC <- rotatePC[which(rotatePC)]
      if(length(rotatePC) >0) pca$x[,rotatePC] <- -1*pca$x[,rotatePC] }
    if(debug) {message(fxNa,"plotPCAw_5 .. lab123: ",wrMisc::pasteC(lab123)); plotPCAw_5 <- list(dat=dat,sampleGrp=sampleGrp,pca=pca,useTit=useTit,cexTxt=cexTxt ) }
    ## start plot
    ch1 <- try(graphics::plot(pca$x[,c(1,2)], main=useTit, cex.axis=0.7*cexTxt, cex.lab=cexTxt*0.75, cex.main=if(ncol(dat) > 3) 1.4 else 0.85, xlab=lab123[1], ylab=lab123[2], las=1, type="n"))   # empty plot
    if(inherits(ch1, "try-error")) {dat <- NULL; message(fxNa,"UNABLE TO DRAW PLOT !!  check current plotting device...")} }

  if(length(dat) >0) {
    ## check for legend-location (if legend should be drawn)
    if(is.character(showLegend) & length(showLegend)==1) {
      chL <- showLegend %in% c("bottomleft","bottomright","topright","topleft")
      if(!chL) showLegend <- TRUE
    }
    displLeg <- if(is.logical(showLegend)) {if(showLegend) checkForLegLoc(matr=pca$x, sampleGrp=sampleGrp, showLegend=showLegend, suplSpace=5.5, silent=silent,callFrom=callFrom) else list(FALSE)
    } else { if(is.character(showLegend)) list(TRUE, showLegend) else list(FALSE)}
    if(debug) {message(fxNa,"plotPCAw_5b .. displLeg: ",displLeg); plotPCAw_5 <- list(dat=dat,sampleGrp=sampleGrp,pca=pca,useTit=useTit,cexTxt=cexTxt,lab123=lab123 ) }

    if(displLeg[[1]] & showLegend) graphics::legend(displLeg[[2]], pch=useSymb.ori, col=colBase,
      paste(1:length(levels(sampleGrp)),"..",substr(unique(averNa3),1,25)), text.col=colBase,
      cex=cexTxt*max(0.4, round(0.75 -((length(unique(sampleGrp)) %/% 5)/31),3)), xjust=0.5, yjust=0.5)
    graphics::points(pca$x[,c(1,2)], pch=useSymb,col=useCol,cex=0.8)
    if(debug) { message(fxNa,"plotPCAw_6 .. displLeg: ",unlist(displLeg[1:2])) }
    ## prepare for labels on points
    if(length(pointLabelPar) >0) {
      chLe <- sapply(pointLabelPar,length) %in% c(1,ncol(dat))
      if(any(!chLe, na.rm=TRUE) & !silent) message(fxNa,"Some elements of argument 'pointLabelPar' may have odd lengths, they might be discarded")
      if(identical(pointLabelPar,TRUE)) .addTextToPoints(x=pca$x, cex=0.6*cexTxt, adj="auto") else {
        if(length(pointLabelPar) >0) .addTextToPoints(x=pca$x, paramLst=pointLabelPar, cex=0.6*cexTxt, adj="auto") }}
    if(debug) {message(fxNa,"plotPCAw_7 .. displBagPl=",displBagPl); plotPCAw_7 <- list(dat=dat,sampleGrp=sampleGrp,pca=pca,useTit=useTit,cexTxt=cexTxt,lab123=lab123,showLegend=showLegend,displLeg=displLeg,
       nGrp=nGrp,averNa3=averNa3,pointLabelPar=pointLabelPar,colBase=colBase,
      outCoef=outCoef,useSymb2=useSymb2,nGrpForMedian=nGrpForMedian,nGrpLim=nGrpLim,displBagPl=displBagPl,pca.grpMed=pca.grpMed,figNumOffs=figNumOffs, callFrom=callFrom,fxNa=fxNa, silent=silent, debug=debug)}

    ## prepare for bagplot
    for(i in 1:length(levels(sampleGrp))) {
      useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256),
        alpha=if(max(table(sampleGrp)) > 2) signif(0.04 +0.25/nGrp, digits=3) else 0.1)
      dispOut <- length(unlist(pointLabelPar)) <1
      outli <- if(displBagPl) addBagPlot(pca$x[which(as.numeric(wrMisc::naOmit(sampleGrp)) %in% i), 1:2], bagCol=useCol2, outCoef=outCoef,
        ctrPch=useSymb2, returnOutL=TRUE, outlCol=useCol2, addSubTi=length(pointLabelPar >0) , callFrom=fxNa, silent=silent, debug=debug) else NULL
      if(length(outli) >0) outL[[i]] <- outli }
       if(debug) {message(fxNa,"plotPCAw_7.",i," .. displBagPl=",displBagPl) }
    if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {              # display of group-center (number & )
      cexCtr <- if(cexTxt <=1) useCex  +0.2 else useCex
      if(!displBagPl & any(duplicated(sampleGrp), na.rm=TRUE)) {  # add group-centers for those where multiple replicates/points were drawn (if not drawn by addBagPot)
        ## which points have replicates ? set point not needed to draw as pch as NA
        pch1 <- useSymb; col1 <- useCol             # indiv points
        if(!all(duplicated(sampleGrp))) pch1 <- rep(NA,length(pch1))
        if(any(duplicated(sampleGrp), na.rm=TRUE))  pch1[ which(!duplicated(sampleGrp, fromLast=TRUE) & !duplicated(sampleGrp, fromLast=FALSE)) ]  <- NA
        graphics::points(pca.grpMed[,c(1,2)], col=useCol, pch=pch1, cex=1.7)    # draw group-means as same color, but bigger symbol
      }
      graphics::text(pca.grpMed[,c(1,2)] + 1.1*figNumOffs(pca, cols=c(1,2)),                     # group-center names
        as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=cexCtr, font=2,col=1) }
    if(debug) {message(fxNa,"plotPCAw_8") }

    ## subtitle
    graphics::mtext(paste0("n=",nrow(dat)," ",rowTyName, " ;  Samples shown as open symbols",
      if(isTRUE(pointLabelPar)) ", color stars and black numbers for group-centers",
      if(displBagPl) ", groups highlighted as bagplot"), cex=cexSub, line=0.25)
    ##  add more plots to include 3rd PC
    if(suplFig & ncol(dat) > 3) {
      graphics::plot(pca$x[,c(2,3)], pch=useSymb, main="PCA : 2nd and 3rd Component", col=useCol, cex=0.6, cex.axis=0.6*cexTxt, cex.lab=cexTxt*0.7, xlab=lab123[2], ylab=lab123[3], las=1 )
      for(i in 1:length(unique(sampleGrp))) {
        useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256), alpha=if(max(table(sampleGrp)) > 2) signif(0.04 +0.25/nGrp,digits=3) else 0.1)
        if(displBagPl) addBagPlot(pca$x[as.numeric(averNa3) %in% i,c(2,3)], outCoef=outCoef, bagCol=useCol2, ctrPch=useSymb2, returnOutL=FALSE, silent=silent, debug=debug) }
      if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {
        graphics::points(pca.grpMed[,c(2,3)], col=colBase[order(as.numeric(unique(averNa3)))], pch=useSymb2, cex=1.1)
        cexCtr2 <- if(cexTxt <=1) useCex -0.1 else useCex
        graphics::text(pca.grpMed[,c(2,3)]+ 1.7*figNumOffs(pca, cols=c(2,3)), as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=cexCtr2, col=1) }
      graphics::plot(pca$x[,c(1,3)], pch=useSymb, main="PCA : 1st and 3rd Component", col=useCol, cex=0.6, cex.axis=0.6*cexTxt, cex.lab=cexTxt*0.7, xlab=lab123[1], ylab=lab123[3],las=1 )
      for(i in 1:length(unique(sampleGrp))) {
        useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256), alpha=if(max(table(sampleGrp)) > 2) signif(0.04+0.25/nGrp,digits=3) else 0.1)
        if(displBagPl) addBagPlot(pca$x[as.numeric(averNa3) %in% i,c(1,3)], outCoef=outCoef, bagCol=useCol2,ctrPch=useSymb2,returnOutL=FALSE,silent=TRUE) }
      if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {
        graphics::points(pca.grpMed[,c(1,3)],col=colBase[order(as.numeric(unique(averNa3)))], pch=useSymb2, cex=1.1)
        graphics::text(pca.grpMed[,c(1,3)] +1.7*figNumOffs(pca, cols=c(1,3)), as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=useCex-0.1, col=1) }
    }

    if(suplFig) { graphics::plot(pca, main="Screeplot on Variance Captured by the Principal Components", cex.main=if(ncol(dat) > 3) 1.3 else 0.8, cex.axis=0.6*cexTxt )
      graphics::mtext(at=(1.2*(1:length(pca$sdev)) -0.5)[1:min(10,ncol(dat))], as.character(1:length(pca$sdev))[1:min(10,ncol(dat))], cex=0.7, side=1) }
    if(getOutL & !is.null(rownames(pca$x))) return(outL)
  } else warning(fxNa,"OMIT plot !  Unable to calculate principal components !")
}



#' @export
.addTextToPoints <- function(x,paramLst=NULL,labels=NULL,cex=NULL,col=NULL, adj="auto",callFrom=NULL,silent=FALSE){
  ##  add text at given position(s) to (scatter) plot (left top display), eg names describing points
  ## 'x' should be matrix or data.frame with numeric coordinates for plotting labels,
  ## paramLst may be list with additional parameters (priority over separate cex or col), otherwise names displayed will be taken from 'labels' or rownames of coordinates 'x'
  ## 'labels' .. for text to be displayed (has not priority over 'paramLst'), if nothing valid found numbers will be displayed
  ## 'adj' .. values other than 0,0.5 or 1 will lead to 'auto' where text is displayed only to left of sufficient space available
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".addTextToPoints")
  textCex <- textLabel <- textAdj <- textCol <- textOffSet <- NULL
  if(length(dim(x)) !=2) message(fxNa," Trouble ahead, 'x' shound be matrix or data.frame of min 2 columns")
  if(length(adj) >0) textAdj <- paramLst$textAdj <- adj
  if(length(cex) >0) textCex <- cex
  if(length(col) ==nrow(x) | length(col)==1) textCol <- col
  textLabel <- if(length(labels) != length(x) ) names(x) else labels
  if(is.list(paramLst) & length(paramLst) >0) {
    textLabel <- if("labels" %in% names(paramLst)) paramLst$labels else {if(is.null(rownames(x))) 1:nrow(x) else rownames(x)}
    if("textLabel" %in% names(paramLst)) textLabel <- paramLst$textLabel
    if("col" %in% names(paramLst)) textCol <- paramLst$col
    if("textCol" %in% names(paramLst)) textCol <- paramLst$textCol
    if("cex" %in% names(paramLst)) textCex <- paramLst$cex
    if("textCex" %in% names(paramLst)) textCex <- paramLst$textCex
    if("offSet" %in% names(paramLst)) textOffSet <- paramLst$offSet
    if("textOffSet" %in% names(paramLst)) textOffSet <- paramLst$textOffSet
    if("adj" %in% names(paramLst)) textAdj <- paramLst$adj
    if(is.null(names(paramLst)) & is.character(paramLst[[1]])) textLabel <- paramLst[[1]]
    if("textAdj" %in% names(paramLst)) {
      if(!any(c(0,0.5,1) %in% paramLst$textAdj, na.rm=TRUE)) {
        strW <- graphics::strwidth(textLabel, cex=if("cex" %in% names(paramLst)) paramLst$cex else 1) + 0.2
        chLe <- graphics::par("usr")[1] < x[,1] -strW
        paramLst$textAdj <- 0 + chLe }
      textAdj <- paramLst$textAdj }
  } else {if(is.character(paramLst) & length(paramLst)==nrow(x)) textLabel <- paramLst}
  if(length(textLabel) != nrow(x)) textLabel <- rownames(x)
  if(is.null(textLabel)) textLabel <- 1:nrow(x)
  if(length(textCol) <1) textCol <- rep(grDevices::grey(0.8), nrow(x)) else if(length(textCol) ==1) textCol <- rep(textCol, nrow(x))
  if(length(textAdj) <1) textAdj <- 1:nrow(x)
  if(length(textOffSet) <1) {ra <- graphics::par("usr"); textOffSet <- c(2,-4)*signif((ra[c(2,4)] - ra[c(1,3)])/300,3) }
  textOffSet <- matrix(rep(textOffSet, each=nrow(x)), ncol=2)
  if(0 %in% textAdj) textOffSet[,1] <- textOffSet[,1]*(2*textAdj -1)
  if(length(unique(textAdj)) >1) { for(i in unique(textAdj)) { j <- which(textAdj==i)
    graphics::text(x[j,c(1,2)] -textOffSet[j,], labels=textLabel[j], col=textCol[j], cex=textCex, adj=i)}
  } else graphics::text(x[,c(1,2)] -textOffSet, labels=textLabel, col=textCol, cex=textCex, adj=textAdj) }
  
