#' PCA plot with bag-plot to highlight groups
#'
#' Function to plot \href{https://en.wikipedia.org/wiki/Principal_component_analysis}{principal components analysis (PCA)}, 
#' with options to show center and potential outliers for each of the groups (columns of data).
#' One of the specificities of this implementation is the integration of bag-plots to better visualize different groups of points 
#' (if they can be organized so beforehand as distinct groups) :
#' The main body of data is shown as 'bag-plots' (a bivariate boxplot, see \href{https://en.wikipedia.org/wiki/Bagplot}{Bagplot})
#' with different transparent colors to highlight the core part of different groups (if they contain more than 2 values per group).
#' Furthermore, group centers are shown as average or median (see 'nGrpForMedian') with stars & index-number (if <25 groups).
#' Layout is automatically set to 2 or 4 subplots (if plotting more than 2 principal components makes sense).
#' Note : This function uses for calulating PCA \code{\link[stats]{prcomp}} with default \code{center=TRUE} and \code{scale.=FALSE}, (different to princomp() which standardizes by default).
#' Note: \code{NA}-values cannot (by definition) be processed by PCA - all lines with any non-finite values/content (eg \code{NA}) will be omitted !
#' Note : Package RColorBrewer may be used if avaialble.
#' Finally, note that several other packages dedicated to PCA exist, for example \href{https://CRAN.R-project.org/package=FactoMineR}{FactoMineR} offers 
#'  a very wide spectrum of possibiities, in particular for combined numeric and categorical data. 
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
#' @param displBagPl (logical) if \code{TRUE}, show bagPlot (group-center) if >3 points per group otherwise the average-confidence-interval  
#' @param getOutL (logical) return outlyer samples/values 
#' @param showLegend (logical) toggle to display legend 
#' @param nGrpForMedian (integer) decide if group center should be displayed via its average or median value: If group has less than 'nGrpForMedian' values, the average will be used, otherwise the median; if \code{NULL} no group centers will be displayed
#' @param pointLabelPar (character) define formatting for optional labels next to points in main figure (ie PC1 vs PC2); may be \code{TRUE} or list containing elments 'textLabel','textCol','textCex',
#'  'textOffSet','textAdj' for fine-tuning
#' @param rowTyName (character) for subtitle : specify nature of rows (genes, proteins, probesets,...)
#' @param rotatePC (integer) optional rotation (by -1) for figure of the principal components specified by index
#' @param suplFig (logical) to include plots vs 3rd principal component (PC) and Screeplot  
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @return plot and optional matrix of outlyer-data
#' @seealso (used in this function for the PCA underneith:) \code{\link[stats]{prcomp}}, \code{\link[stats]{princomp}}, the package \href{https://CRAN.R-project.org/package=FactoMineR}{FactoMineR} 
#' @examples
#' set.seed(2019); dat1 <- matrix(round(c(rnorm(1000), runif(1000,-0.9,0.9)),2), 
#'   ncol=20, byrow=TRUE) + matrix(rep(rep(1:5,6:2), each=100), ncol=20)
#' biplot(prcomp(dat1))      # traditional plot
#' (grp = factor(rep(LETTERS[5:1],6:2)))
#' plotPCAw(dat1, grp)
#' @export
plotPCAw <- function(dat, sampleGrp, tit=NULL, useSymb=c(21:25,9:12,3:4), center=TRUE, scale.=TRUE, colBase=NULL, useSymb2=NULL,
  displBagPl=TRUE, getOutL=FALSE, cexTxt=1, showLegend=TRUE, nGrpForMedian=6, pointLabelPar=NULL, rowTyName="genes",
  rotatePC=NULL, suplFig=TRUE, callFrom=NULL, silent=FALSE) {
  ## note : so far not well adopted for plotting all points in transparent grey with same symb for interactive mouse-over
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotPCAw")
  msg <- " 'dat' must be matrix, data.frame (with at least 1 line and 2 columns) or list/MArrayLM-object (with element 'datImp','dat' or 'data' as matrix or data.frame)"
  if("MArrayLM" %in% class(dat) | is.list(dat)) { chDat <- c("datImp","dat","data","quant","datRaw")  %in% names(dat)
    chGr <- c("groups","grp","sampleGroups")  %in% names(dat)
    if(any(chGr) & length(sampleGrp) <1) sampleGrp <- dat[[which(names(dat) %in% c("groups","grp","sampleGroups"))[1]]] 
    if(any(chDat)) dat <- dat[[which(names(dat) %in% c("datImp","dat","data","datRaw"))[1]]] else stop(msg) }
  if(any(length(dim(dat)) !=2, dim(dat) < 1:2)) stop(msg)
  if(length(sampleGrp) <1 | sum(is.na(sampleGrp))==length(sampleGrp)) stop("argument 'sampleGrp' seems to be empty (or all NA)")
  if(length(sampleGrp) != ncol(dat)) stop(fxNa," length of argument 'sampleGrp' must match number of columns of 'dat'")
  if(!is.factor(sampleGrp)) sampleGrp <- as.factor(sampleGrp)
  nGrpLim <- 25                      # decide/limit if group-centers are shown as stars
  ## reduce to columns defined by 'sampleGrp' (groups declared as NA will be removed) and to columns with finite data
  dat <- dat[,!is.na(sampleGrp)]
  chFin <- is.finite(dat)
  chCol <- colSums(chFin) <2
  if(sum(chCol) >0){
    dimIni <- dim(dat)
    if(!silent) message(fxNa,"removing ",sum(chCol)," (out of ",ncol(dat),") columns with less than 2 finite values")
    dat <- dat[,which(!chCol)]
    chFin <- chFin[,which(!chCol)]
    if(ncol(dat) <2) stop("after removing bad (non-finite) columns less than 2 columns remain !")
    sampleGrp <- sampleGrp[which(!chCol)]
    if(length(useSymb) == dimIni[2]) useSymb <- useSymb[which(!chCol)] }
  opar <- graphics::par(no.readonly=TRUE) 
  on.exit(graphics::par(opar))  
  if(suplFig & ncol(dat) > 3) on.exit(graphics::par(opar))  else {
    on.exit(graphics::par(opar$mar))
    on.exit(graphics::par(opar$mfrow))
  }
  ## remove rows with NAs 
  chRow <- rowSums(!chFin)                      # check for lines with non-finite values (avoid error during svd/pca ..)
  if(any(chRow >0)) if(sum(chRow >0) ==nrow(dat)) stop("all lines have some non-finite values !?!") else {
    if(!silent) message(fxNa," eliminating ",sum(chRow >0)," (out of ",nrow(dat),") lines with non-finite values")
    dat <- dat[which(chRow <1),] }
  averNa3 <- wrMisc::naOmit(sampleGrp)
  if(suplFig & ncol(dat) > 3) graphics::layout(matrix(c(1,1:3,4,4),3,2,byrow=TRUE), heights=c(10,5,4)) 
  if(!identical(ncol(dat),length(wrMisc::naOmit(sampleGrp))) & !silent) warning("'plotPCA': number of columns and number of (QC-filtered) arrays in 'AverNa2' don't match !")
  nGrp <- length(levels(sampleGrp))
  if(!is.null(colBase)) if(length(colBase) != nGrp) {
    if(!silent) message(fxNa," ignore invalid 'colBase', expecting length ",nGrp)
    colBase <- NULL }
  chPa <- requireNamespace("RColorBrewer", quietly=TRUE)
  if(!silent & !chPa) message(fxNa,"Please install package 'RColorBrewer' for improved color-gradients") 
  if(is.null(colBase)) colBase <- if(nGrp <10 & chPa) {RColorBrewer::brewer.pal(max(3,nGrp),"Set1")[1:nGrp]
    } else {
    grDevices::rainbow(1.5 +nGrp*1.08)[1:(nGrp+1)][-ceiling(length(levels(sampleGrp))/2)]}
  useCol <- colBase[order(as.numeric(unique(averNa3)))]
  if(!identical(unique(table(averNa3)),as.integer(1))) useCol <- useCol[as.numeric(wrMisc::naOmit(sampleGrp))]
  ## MAIN
  pca <- stats::prcomp(t(dat),center=center,scale.=scale.)
  percVar <- 100*round(pca$sdev^2/sum(pca$sdev^2),3)
  pca.grpMed <- pca.grpAv <- matrix(nrow=nGrp, ncol=min(4,ncol(pca$x)))
  tabAverNa3 <- table(averNa3)
  if(length(nGrpForMedian) >0) { for(i in 1:min(4,ncol(pca$x))) {
      pca.grpMed[,i] <- unlist(tapply(pca$x[,i], averNa3, stats::median, na.rm=TRUE))
      pca.grpAv[,i] <- unlist(tapply(pca$x[,i], averNa3, mean, na.rm=TRUE)) }
    rownames(pca.grpMed) <- rownames(pca.grpAv) <- unique(averNa3)
    pca.grpMed[names(tabAverNa3)[which(tabAverNa3 <nGrpForMedian)],] <- pca.grpAv[names(tabAverNa3)[which(tabAverNa3 <nGrpForMedian)],]}
  useCex <- round( (1/nGrp +0.65)*cexTxt, 2)
  useSymb.ori <- rep(useSymb, 1 +length(levels(sampleGrp)) %/% length(useSymb)) [1:length(levels(sampleGrp))] #
  useSymb <- (useSymb.ori[order(unique(as.numeric(averNa3)))] )[as.numeric(wrMisc::naOmit(sampleGrp))]
  useTxtCol <- 1
  figNumOffs <- function(pca, cols=c(1,2), fact=100) signif(c(abs(diff(range(pca$x[,cols[1]])))/fact,abs(diff(range(pca$x[,cols[2]])))/fact), digits=3)
  outL <- list()
  length(outL) <- length(unique(averNa3)); names(outL) <- unique(averNa3)
  lab123 <- paste("PC",1:3," (",percVar[1:3],"%)",sep="")
  useTit <- if(is.null(tit)) "Principal Components of Samples" else tit
  if(suplFig) useTit <- paste(useTit, if(nchar(useTit) < 35) ": 1st and 2nd Component" else ", PC1 & PC2")
  ## optional rotate axis
  if(length(rotatePC) >0 & is.numeric(rotatePC)) {
    chRo <- rotatePC >0 & rotatePC <= ncol(dat)
    if(any(!rotatePC)) rotatePC <- rotatePC[which(rotatePC)]
    if(length(rotatePC) >0) pca$x[,rotatePC] <- -1*pca$x[,rotatePC]
  }  
  ## start plot
  graphics::plot(pca$x[,c(1,2)], main=useTit, type="n", cex.axis=0.7*cexTxt, cex.lab=cexTxt*0.75, cex.main=if(ncol(dat) > 3) 1.4 else 0.85, xlab=lab123[1], ylab=lab123[2], las=1)
  displLeg <- checkForLegLoc(matr=pca$x, sampleGrp=sampleGrp, showLegend=showLegend, suplSpace=5.5, callFrom=callFrom)
  if(displLeg[[1]]) graphics::legend(displLeg[[2]], pch=useSymb.ori, col=colBase,
    paste(1:length(levels(sampleGrp)),"..",substr(unique(averNa3),1,25)), text.col=colBase,
    cex=cexTxt*max(0.4, round(0.75 -((length(unique(sampleGrp)) %/% 5)/31),3)), xjust=0.5, yjust=0.5)
  graphics::points(pca$x[,c(1,2)], pch=useSymb,col=useCol,cex=0.8)
  ## prepare for labels on points
  if(length(pointLabelPar) >0) {
    chLe <- sapply(pointLabelPar,length) %in% c(1,ncol(dat))
    if(any(!chLe) & !silent) message(fxNa," some elements of argument 'pointLabelPar' may have odd lengths, they might be discarded") 
    if(identical(pointLabelPar,TRUE)) .addTextToPoints(x=pca$x, cex=0.6*cexTxt)
    if(length(pointLabelPar) >0) .addTextToPoints(x=pca$x, paramLst=pointLabelPar, cex=0.6*cexTxt) } 
  ## prepare for bagplot
  for(i in 1:length(levels(sampleGrp))) {                     
    useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256),
      alpha=if(max(table(sampleGrp)) > 2) signif(0.04 +0.25/nGrp, digits=3) else 0.1)
    dispOut <- length(unlist(pointLabelPar)) <1
    outli <- if(displBagPl) addBagPlot(pca$x[as.numeric(wrMisc::naOmit(sampleGrp)) %in% i, 1:2], outCoef=3.5, bagCol=useCol2,
      ctrPch=useSymb2, returnOutL=TRUE, outlCol=useCol2, callFrom=fxNa, silent=TRUE) else NULL
    if(length(outli) >0) outL[[i]] <- outli }
  if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {              # display of group-center (number & )
    cexCtr <- if(cexTxt <=1) useCex+0.2 else useCex
    if(!displBagPl) graphics::points(pca.grpMed[,c(1,2)], col=colBase[order(as.numeric(unique(averNa3)))], pch=useSymb2, cex=1.7)  # drawn by addBagPlot
    graphics::text(pca.grpMed[,c(1,2)] + 1.1*figNumOffs(pca, cols=c(1,2)),                     # group-center names
      as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=cexCtr, font=2,col=1) }
  ## subtitle
  graphics::mtext(paste("n=",nrow(dat)," ",rowTyName," ;  Samples shown as open symbols",
    if(identical(pointLabelPar,TRUE)) ", color stars and black numbers for group-centers",
    if(displBagPl) ", groups highlighted as bagplot", sep=""), cex=0.6, line=0.25)
  ##  add more plots to include 3rd PC
  if(suplFig & ncol(dat) > 3) {
    graphics::plot(pca$x[,c(2,3)], pch=useSymb, main="PCA : 2nd and 3rd Component", col=useCol, cex=0.6, cex.axis=0.6*cexTxt, cex.lab=cexTxt*0.7, xlab=lab123[2], ylab=lab123[3],las=1 )
    for(i in 1:length(unique(sampleGrp))) {
      useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256), alpha=if(max(table(sampleGrp)) > 2) signif(0.04+0.25/nGrp,digits=3) else 0.1)
      if(displBagPl) addBagPlot(pca$x[as.numeric(averNa3) %in% i,c(2,3)], outCoef=3.5, bagCol=useCol2,ctrPch=useSymb2,returnOutL=FALSE,silent=TRUE) }
    if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {
      graphics::points(pca.grpMed[,c(2,3)], col=colBase[order(as.numeric(unique(averNa3)))], pch=useSymb2, cex=1.1)
      cexCtr2 <- if(cexTxt <=1) useCex -0.1 else useCex
      graphics::text(pca.grpMed[,c(2,3)]+ 1.7*figNumOffs(pca, cols=c(2,3)), as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=cexCtr2, col=1) }
    graphics::plot(pca$x[,c(1,3)], pch=useSymb, main="PCA : 1st and 3rd Component", col=useCol, cex=0.6, cex.axis=0.6*cexTxt, cex.lab=cexTxt*0.7,xlab=lab123[1],ylab=lab123[3],las=1 )
    for(i in 1:length(unique(sampleGrp))) {
      useCol2 <- grDevices::rgb(t(grDevices::col2rgb(colBase[order(as.numeric(unique(averNa3)))][i])/256), alpha=if(max(table(sampleGrp)) > 2) signif(0.04+0.25/nGrp,digits=3) else 0.1)
      if(displBagPl) addBagPlot(pca$x[as.numeric(averNa3) %in% i,c(1,3)], outCoef=3.5, bagCol=useCol2,ctrPch=useSymb2,returnOutL=FALSE,silent=TRUE) }
    if(length(nGrpForMedian) >0 & length(levels(sampleGrp)) < nGrpLim & !is.null(useSymb2)) {
      graphics::points(pca.grpMed[,c(1,3)],col=colBase[order(as.numeric(unique(averNa3)))], pch=useSymb2, cex=1.1)
      graphics::text(pca.grpMed[,c(1,3)] +1.7*figNumOffs(pca, cols=c(1,3)), as.character(1:nGrp)[order(as.numeric(unique(averNa3)))], cex=useCex-0.1, col=1) }
  }
  if(suplFig) { graphics::plot(pca, main="Screeplot on Variance Captured by the Principal Components", cex.main=if(ncol(dat) > 3) 1.3 else 0.8, cex.axis=0.6*cexTxt )
    graphics::mtext(at=(1.2*(1:length(pca$sdev)) -0.5)[1:min(10,ncol(dat))], as.character(1:length(pca$sdev))[1:min(10,ncol(dat))], cex=0.7, side=1) }
  if(getOutL & !is.null(rownames(pca$x))) return(outL)
  }

#' @export
.addTextToPoints <- function(x,paramLst=NULL,labels=NULL,cex=NULL,col=NULL,callFrom=NULL,silent=FALSE){
  ##  add text at given position(s) to (scatter) plot (left top display), eg names describing points
  ## 'x' should be matrix or data.frame with numeric coordinates for plotting labels, 
  ## paramLst may be list with additional parameters (priority over separate cex or col), otherwise names displayed will be taken from 'labels' or rownames of coordinates 'x'
  ## 'labels' .. for text to be displayed (has not priority over 'paramLst'), if nothing valid found numbers will be displayed
  fxNa <- wrMisc::.composeCallName(callFrom,newNa=".addTextToPoints")
  textCex <- textLabel <- textAdj <- textCol <- textOffSet <- NULL
  if(length(cex) >0) textCex <- cex
  if(length(col) ==nrow(x) | length(col)==1) textCol <- col
  if(length(labels) ==nrow(x)) textLabel <- labels
  if(is.list(paramLst)) {
    if("labels" %in% names(paramLst)) textLabel <- paramLst$labels
    if("textLabel" %in% names(paramLst)) textLabel <- paramLst$textLabel
    if("col" %in% names(paramLst)) textCol <- paramLst$col
    if("textCol" %in% names(paramLst)) textCol <- paramLst$textCol
    if("cex" %in% names(paramLst)) textCex <- paramLst$cex
    if("textCex" %in% names(paramLst)) textCex <- paramLst$textCex
    if("offSet" %in% names(paramLst)) textOffSet <- paramLst$offSet
    if("textOffSet" %in% names(paramLst)) textOffSet <- paramLst$textOffSet
    if("textAdj" %in% names(paramLst)) textAdj <- paramLst$textAdj
    if("adj" %in% names(paramLst)) textAdj <- paramLst$adj
    if(is.null(names(paramLst)) & is.character(paramLst[[1]])) textLabel <- paramLst[[1]]  
  } else {if(is.character(paramLst) & length(paramLst)==nrow(x)) textLabel <- paramLst}
  if(length(textLabel) != nrow(x)) textLabel <- rownames(x)
  if(is.null(textLabel)) textLabel <- 1:nrow(x)
  if(length(textCol) <1) textCol <- grDevices::grey(0.85)
  if(length(textAdj) <1) textAdj <- 1:nrow(x)
  if(length(textOffSet) <1) textOffSet <- c(1,-4)*signif(apply(x[,1:2], 2, function(x) abs(diff(range(x))))/300,3)   # design for left top display
  textOffSet <- matrix(rep(textOffSet, each=nrow(x)), ncol=2)
  graphics::text(x[,c(1,2)]-textOffSet, labels=textLabel, col=textCol, cex=textCex, adj=textAdj) }
   
