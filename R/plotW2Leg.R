#' x-y plot with 2 legends  
#'
#' This is a modified version of \code{\link[base]{plot}} for 2-dimensional data,
#' allowing to choose symbols and colors of points according to an additional columns of \code{dat}. 
#'   
#' 
#'
#' @param dat (matrix or data.frame) main input
#' @param useCol (character or integer) columns form \code{dat}: The 1st and 2nd column are used as x- and y-axis
#' @param tit (character) optional custom title
#' @param subTi (character) optional custom subtitle
#' @param subCex (numeric) cex-like expansion factor for subtitle (see also \code{\link[graphics]{par}})
#' @param pch (integer) symbols to use for plotting (see also \code{\link[graphics]{par}}), will be associated to 4th column of \code{useCol}
#' @param xlim (numeric, length=2) x- axis limits (see also \code{\link[graphics]{par}})
#' @param ylim (numeric, length=2) y- axis limits (see also \code{\link[graphics]{par}})
#' @param xlab (character) custom x-axis label
#' @param ylab (character) custom x-axis label 
#' @param ablines (list) optional horzontal and/or vertical gray dashed guide-lines
#' @param legendloc (character) location of legend (of symbols) 
#' @param txtLegend (character) optional label for legend (of symbols) 
#' @param histLoc (character) location of histomgram-legend (of 3rd column of \code{useCol}) 
#' @param legHiTi (character) optional title for histomgram-legend
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @seealso \code{\link[base]{plot}}   
#' @return graphical output only
#' @examples
#' x1 <- cbind(x=c(2,1:7), y=8:1 +runif(8), grade=rep(1:4,2))
#' plotW2Leg(x1,useCol=c("x","y","y","grade"))
#' @export
plotW2Leg <- function(dat,useCol=c("logp","slope","medAbund","startFr"), tit=NULL, subTi=NULL, subCex=0.9, pch=21:25, xlim=NULL, ylim=NULL,
  xlab=NULL, ylab=NULL, ablines=NULL, legendloc="topright", txtLegend=NULL, histLoc="bottomleft", legHiTi=NULL, silent=TRUE, callFrom=NULL) {
  ## plot with 2 extra legends
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotW2Leg")
  if(is.null(tit)) tit <- "log p-value vs slope"
  msg <- "Invalid entries for 'useCol' (should be column-names or integer as index)"
  if(length(useCol) <4) useCol <- c(useCol, rep(NA,4 -length(useCol)))
  if(!is.integer(useCol)) {chCol <- match(useCol, colnames(dat))
    if(all(is.na(chCol))) stop(msg)
    if(any(is.na(chCol)) & !silent) message(fxNa,"Trouble ahead ! Can't find columns ",wrMisc::pasteC(useCol[is.na(chCol)],quoteC="'"))
    useCol <- chCol 
  } else if(any(useCol <1 | useCol >ncol(dat))) stop(msg)
  useColor <- if(!is.na(useCol[3])) wrMisc::colorAccording2(dat[,useCol[3]], gradTy="rainbow", revCol=TRUE, nEndOmit=14) else grDevices::gray(0.6)
  stRa <- if(!is.na(useCol[4])) range(dat[,useCol[4]], na.rm=TRUE) else NA
  graphics::plot(dat[,useCol[1:2]], type="n", main=tit, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, las=1) 
  if(length(subTi) >0) graphics::mtext(subTi, cex=subCex)
  ## guide-lines
  if(length(ablines) >0) { 
    if("v" %in% names(ablines)) graphics::abline(v=ablines$v, lty=2, col="grey") 
    if("v" %in% names(ablines)) graphics::abline(h=ablines$v, lty=2, col="grey") }
  ## legend for symbols
  if(!is.na(useCol[4]) & length(legendloc)==1) graphics::legend(legendloc, paste(                                                                  
    if(length(txtLegend) <1) paste(colnames(dat)[useCol[4]],"=") else txtLegend,
    stRa[1]:stRa[2]), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
  if(length(subTi) >0) graphics::mtext(subTi,cex=0.9)
  ## histigram legend 
  if(!is.na(useCol[3]) & length(histLoc)==1) {
    if(length(legHiTi) <1) legHiTi <- colnames(dat)[useCol[3]]
    hi1 <- graphics::hist(dat[,useCol[3]], plot=FALSE)
    legendHist(sort(dat[,useCol[3]]), colRamp=useColor[order(dat[,useCol[3]])][cumsum(hi1$counts)], location=histLoc, legTit=legHiTi) }
  ## the main points (on top)
  pch1 <- if(!is.na(useCol[4])) pch[dat[,useCol[4]]] else {if(length(pch)==nrow(dat)) pch else rep(pch[1], nrow(dat))}
  ptCol <- useColor
  if(any(pch1 %in% 21:25)) ptCol[which(pch1 %in% 21:25)] <- 1 
  graphics::points(dat[,useCol[1:2]], col=ptCol, bg=useColor, pch=pch1) 
  }
  
