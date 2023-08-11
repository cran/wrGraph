#' Plot profiles according to CLustering
#'
#' This function was made for visualuzing the result of clustering of a numeric vector or clustering along multiple columns of a matrix.
#' The data will be plotted like a reglar scatter-plot, but some extra space is added to separate clusters and dashed lines highlight cluster-borders.
#' If no mean/representative value is spacified, a geometric mean will be calculated along all columns of  \code{dat}.
#' In case \code{dat} has multiple columns, a legend and a representative (default geometric mean) dashed grey line will be displayed.
#' 
#' @param dat (matrix or data.frame) main input with data to plot as points
#' @param clu (numeric or character) clustering results; if length=1 and character this term will be understood as colum-name with cluster-numbers from \code{dat}  
#' @param meanD (numeric) mean/representative of multiple series for display as lines; if length=1 and character this term will be understood as columname with cluster-numbers from \code{dat}
#' @param tit (character) optional custom title
#' @param col (character) custom colors
#' @param pch (integer) custom plotting symbols (see also \code{\link[graphics]{par}})
#' @param xlab (character) custom x-axis label
#' @param ylab (character) custom y-axis label
#' @param cex (numeric) cex-like expansion factor  (see also \code{\link[graphics]{par}})
#' @param cexTit (numeric) cex-like expansion factor for title  (see also \code{\link[graphics]{par}})

#' @param meCol (character) color for (dashed) line of mean/representative values
#' @param meLty (integer) line-type line of mean/representative values (see also \code{lty} in \code{\link[graphics]{par}})
#' @param meLwd (numeric) line-width line of mean/representative values (see also \code{lwd} in \code{\link[graphics]{par}})
#' @param legLoc (character) legend location
#' @param silent (logical) suppress messages
#' @param debug (logical) additonal messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This functin returns a plot only
#' @examples
#' set.seed(2020); dat1 <- runif(12)/2 + rep(6:8, each=4)
#' dat1Cl <- stats::kmeans(dat1, 3)$cluster
#' dat1Cl <- 5- dat1Cl              # bring cluster-numbers in ascending form
#' dat1Cl[which(dat1Cl >3)] <- 1    # bring cluster-numbers in ascending form
#' profileAsClu(dat1, clu=dat1Cl)
#' @export 
profileAsClu <- function(dat, clu, meanD=NULL, tit=NULL, col=NULL, pch=NULL, xlab=NULL, ylab=NULL, meCol="grey", meLty=1, meLwd=1, cex=NULL, cexTit=NULL, legLoc="bottomleft", silent=TRUE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="profileAsClu")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE

  argNames <- c(deparse(substitute(dat)), deparse(substitute(clu)), deparse(substitute(z)))
  msg <- "Invalid argument 'dat'; must be matrix or data.frame with min 2 lines"
  if(length(dat) <1) stop(msg)
  if(length(dim(dat)) <2) dat <- as.matrix(dat)
  if(any(length(dim(dat)) !=2, dim(dat) <1)) stop(msg)
  if(length(clu)==1) { if(is.character(clu)) cluC <- which(colnames(dat)==clu)
    clu <- dat[,cluC]; dat <- dat[,-cluC] }
  if(length(clu) != nrow(dat)) stop("Length of 'clu' does not match number of lines from 'dat')")  
  chClu <- clu[-1] - clu[-length(clu)]
  if(any(any(chClu <0) && any(chClu >0), any(chClu < -1), any(chClu >1))) stop("Invalid 'clu'; cluster numbers must be sorted") 
  if(debug) message(fxNa,"pAC1")
  ## main
  nClu <- length(unique(clu))
  ## geometrix mean
  if(length(meanD)==1) { if(is.character(meanD)) meanC <- which(colnames(dat)==meanD)
    if(length(meanC)==1) {meanD <- dat[,meanC]; dat <- dat[,-meanC]} else {
      message(fxNa,"Can't find column '",meanD,"'")
      meanD <- NULL} }
  if(length(dim(dat)) <2) dat <- as.matrix(dat)
  if(is.null(meanD) && ncol(dat) >1) meanD <- apply(dat, 1, prod)^(1/ncol(dat))
  cluLoc <- table(clu)[rank(unique(clu))]
  cluBord <- cumsum(cluLoc[-nClu]) +0.5
  cluLoc <- cumsum(cluLoc) - cluLoc/2 +0.5
  meanLi <- matrix(NA, nrow=nrow(dat) +nClu, ncol=2)
  inc <- 0
  if(debug) message(fxNa,"pAC2")
  for(i in unique(clu)) { li <- which(clu==i); meanLi[inc+li,] <- cbind(li, meanD[li]); inc <- inc +1}
  meanLi[which(is.na(meanLi[,1])),1] <- c(cluBord,NA)
  if(is.null(col)) col <- 1:ncol(dat)
  if(is.null(pch)) pch <- unique(clu)
  if(is.null(ylab)) ylab <- argNames[1]
  ch1 <- try(graphics::plot(c(range(dat,na.rm=TRUE), rep(NA,nrow(dat)-2)), type="n", las=1, ylab=ylab, main=tit, cex=cex, cex.main=cexTit, xlab=xlab), silent=TRUE)
  if(inherits(ch1, "try-error")) warning(fxNa,"UNABLE TO PRODUCE PLOT !  \n  initial message : ",ch1) else {
    graphics::mtext(paste("clu", unique(clu)), side=3, at=cluLoc)                    # cluster names
    if(ncol(dat) >1) graphics::lines(meanLi, col=meCol, lty=meLty, lwd=meLwd)        # geom mean line
    graphics::abline(v=cluBord, lty=4, col=grDevices::grey(0.8))                     # clu borders
    for(i in 1:ncol(dat)) graphics::points(1:nrow(dat), dat[,i], pch=pch[clu], col=col[i])
    if(ncol(dat) >1) ch1 <- try(graphics::legend(legLoc,c(colnames(dat),"geomMean"), text.col=c(col,1), pch=c(rep(1,ncol(dat)), NA),
      lty=c(rep(NA,ncol(dat)),meLty), lwd=c(rep(NA,ncol(dat)),meLwd), col=c(col[1:ncol(dat)], meCol), cex=0.85, xjust=0.5, yjust=0.5), silent=TRUE)
	 if(inherits(ch1, "try-error")) message(fxNa,"Unable to plot legend")  }
}
   
