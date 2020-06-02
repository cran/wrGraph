#' Separate and plot data by 2 groups 
#'
#' Plot series of data as membership of 2 different grouping vectors (eg by grp=patient and grp2=age-group). 
#'    
#' @param dat (numeric) main data (may contain \code{NA})
#' @param grp (character or factor) grouping of columns of 'dat', eg replicate association
#' @param grp2 (character or factor) aadditional/secondary grouping of columns of 'dat'
#' @param col (character or integer) use custom colors, see also \code{\link[graphics]{par}} 
#' @param pch (integer) symbol to mark group-center  (see also \code{\link[graphics]{par}}) 
#' @param tit (character) custom title
#' @param cex (numeric) expansion factor for text (see also \code{\link[graphics]{par}}) 
#' @param lwd (integer) line-width  (see also \code{\link[graphics]{par}}) 
#' @param lty (integer) line-type  (see also \code{\link[graphics]{par}}) 
#' @param yLab (character) custom y-axis label 
#' @param cexLab (numeric) expansion factor for labels: 1st value for main groups (\code{grp}, eg genotypes), 2nd for detailed text (\code{grp2}, eg animal IDs) (see also \code{\link[graphics]{par}}) 
#' @param sepLines (logical) optional drawing of horizontal lines aiming to separate groups (in analogy to support vectors)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations, or returns data.frame with annot and quant if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) 
#' @examples
#' set.seed(2020); rand1 <- round(runif(12),2) +rep(1:3,each=4)
#' plotBy2Groups(rand1,gl(2,6,labels=LETTERS[5:6]),gl(4,3,labels=letters[1:4]))
#'  
#' @export
plotBy2Groups <- function(dat,grp,grp2=NULL,col=NULL,pch=NULL,tit=NULL,cex=2,lwd=0.5,lty=2,yLab=NULL,cexLab=NULL,sepLines=FALSE,silent=FALSE,callFrom=NULL) {
  ## plot indiv values as membership of 2 grouping vectors (eg by grp=patient and grp2=age-group)
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="plotBy2Groups")  
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  namesXYZ <- c(deparse(substitute(dat)),deparse(substitute(grp)),deparse(substitute(grp2)))
  ## check for undefined elments
  if(length(grp2) <1) grp2 <- grp
  dat0 <- dat <- as.numeric(dat)
  chNaGrp <- is.na(grp) | is.na(grp2)
  if(any(chNaGrp)) { if(all(chNaGrp)) stop(" No individual defined by both groups")
    if(!silent) message(fxNa," remove ",sum(chNaGrp)," due to NAs in (one of the) groups")
    dat <- dat[which(!chNaGrp)]
    dat0 <- dat0[which(!chNaGrp)]
    grp <- grp[which(!chNaGrp)]
    grp2 <- grp2[which(!chNaGrp)]
  }
  ## add'l step to check if grp2 is in correct order ... ?
  zz <- as.numeric(as.factor(grp))
  chBr <- sum(zz[-1] - zz[-length(zz)] !=0) >= length(unique(zz))
  if(chBr) {   # update order of grp
    chOrd1 <- order(grp)
    grp <- grp[chOrd1]
    grp2 <- grp2[chOrd1]
    dat <- dat[chOrd1] 
    dat0 <- dat0[chOrd1] 
  }
  if(length(grp2) < length(dat0)) stop(" 'grp2' not matching 'dat'")
  ## work with grp2 (eg group of patients)
  nGrp2 <- table(grp2)[rank(unique(grp2))]
  grp2ctr <- cumsum(nGrp2) -sapply(nGrp2,function(x) mean(0:(x-1),na.rm=TRUE))
  grp2b <- rep(1:length(unique(grp2)),nGrp2)                  # same regrouping but starts at 1 and increases by 1 
  ## organize dat by grp2   
  dat <- by(as.numeric(dat),grp2b,as.numeric)             # ranking remaines
  names(dat) <- names(nGrp2)
  ## work with grp (eg group 'age-class', already organized by grp2) 
  grp3 <-  tapply(grp2,grp,function(x) unique(x))
  grp3 <- grp3[match(unique(grp),names(grp3))]
  nGrp3 <- sapply(grp3,length)
  grp3ctr <- cumsum(nGrp3) -sapply(nGrp3,function(x) mean(0:(x-1),na.rm=TRUE))
  ## prepare for plot
  if(is.null(col)) col <- as.numeric(as.factor(grp))
  if(is.null(pch)) pch <- as.numeric(as.factor(grp))
  if(length(cexLab) <2) cexLab <- c(0.9,0.7)
  if(is.null(tit)) tit <- paste(namesXYZ[1],"organized by",if(sum(nchar(namesXYZ[2:3])) >19) "two factors" else paste(namesXYZ[2],"and",namesXYZ[3]))
  graphics::plot(grp2b,unlist(dat),col=col,pch=pch,main=tit,cex=2,las=1,xaxt='n',xlab="",ylab=yLab)
  graphics::mtext(at=unique(grp2b),names(grp2ctr),cex=0.7,side=1,col=col[grp2ctr])      # animal no
  graphics::mtext(at=grp3ctr,names(grp3ctr),cex=0.9,side=1,line=1.5,col=unique(col))    # genotype
  ## separation lines
  if(sepLines) {
    datG <- by(dat0,grp,as.numeric)
    datG <- datG[match(unique(grp),names(datG))]
    datGr <- sapply(datG,range,na.rm=TRUE)
    ra0 <- datGr[,-1] - datGr[2:1,-ncol(datGr)]
    if(length(dim(ra0)) <2) ra0 <- as.matrix(ra0)
    ch <- if(length(dim(ra0)) >1) abs(ra0[1,]) < abs(ra0[2,]) else abs(ra0[1]) < abs(ra0[2]) 
    sepLi <- apply(cbind(1:(ncol(datGr)-1),1+ch),1,function(x) datGr[x[2],x[1]] + ra0[3-x[2],x[1]]/2)
    graphics::abline(h=sepLi,col=unique(col)[1:length(sepLi)],lty=lty,lwd=lwd)
    }
}
  
