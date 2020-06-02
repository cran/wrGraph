#' Staggered Chart for Ploting Counts to Multiple Leveles of the Threshold used 
#'
#' The basic idea of this plot is to show how counts data change while shifting a threshold-criterium.
#' At each given threshold the counts are plotted like a staggered bar-chart (or staggered histogram) but without vertical lines to illustrated the almost continuous change
#' from preceedig or following threshold-value.
#' Initially this plot was designed for showing the absolute count-data used when constructing roc-curves (eg using 
#' the package package \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo} with \code{\link[wrProteo]{summarizeForROC}} ).
#' The main input should furnish the panel of threshold as one column and the coresponding counts data as min 2 columns. 
#' The threshold coumns gets specified using the argument \code{threColumn}, the counts-data may either be specified using argument \code{countsCol}
#' or be searched using \code{\link[base]{grep}} using column-names containing the text given in argument \code{varCountNa} with may be combined with 
#' a fixed preceeding part given as argument \code{fixedCountPat}.
#'
#' Investigate count data prepared for plotting ROC curves : cumulative counts plot by species (along different statistical test thresholds).
#' Note : Package \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo} may be used to prepare input (matrix of ROC data).
#'	 
#' @param roc (numeric matrix or data.frame) main input: one column with thresholds and multiple columns of assoicated count data    
#' @param threColumn (integer or character) to specify the column with threshold-data, in typica proteomics benchmark studies this would be 'alph' (for the statistical test threshold)
#' @param countsCol (character of integer, min length=2) choice of column(s) with count-data in 'roc' to be used for display, if not \code{NULL} will override alternative search of columns using 'varCountNa' and 'fixedCountPat'
#' @param varCountNa (character) alternative way to select the columns from 'roc': searched using \code{\link[base]{grep}} using column-names containing the text given in argument \code{varCountNa} with may be combined with a fixed preceeding part given as argument \code{fixedCountPat}
#'  In proteomics benchmark studies this would typically be the species-abbreciations (eg 'H','S','E')
#' @param fixedCountPat (character) optional pattern to help identifying counts-data: if not \code{NULL} it will be used as fixed part in column names to get pasted to \code{varCountNa}.
#'  In proteomics benchmark studies this would typically be 'n.pos.'
#' @param sortAscending (logical) decide if data should be sorted ascending or descending
#' @param vertLine (numeric) for optional vertical line, typically used to highlight alpha 0.05
#' @param col (character) custom colors, see also \code{\link[graphics]{par}}
#' @param tit (character) cutom title
#' @param logScale (logical) display threshld values (x-axis) on log-scale
#' @param las.alph (numeric) orientation of label of alpha-cutoff, see also \code{\link[graphics]{par}}
#' @param displMaxSpec (logical) display on right side of figure max count value of contributing group species 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return plot only
#' @seealso \code{\link[stats]{ecdf}}, for preparing input to ROC: function \code{\link[wrProteo]{summarizeForROC}} in package \href{https://CRAN.R-project.org/package=wrProteo}{wrProteo}
#' @examples
#' set.seed(2019); test1 <- cbind(a=sample.int(n=7,size=50,repl=TRUE),
#'   b=sample.int(n=11,size=50,repl=TRUE),c=sample.int(n=18,size=50,repl=TRUE))
#' test1 <- cbind(alph=seq(0,1,length.out=50),a=cumsum(test1[,1]),b=cumsum(test1[,2]),
#'   c=cumsum(test1[,3]))
#' staggerdCountsPlot(test1,countsCol=c("a","b","c"))
#' ## example below requires the package wrProteo
#' @export
staggerdCountsPlot <- function(roc, threColumn=1, countsCol=NULL, fixedCountPat="n.pos.", varCountNa=NULL, sortAscending=TRUE, vertLine=NULL,
  col=NULL, tit=NULL, logScale=FALSE, las.alph=2, displMaxSpec=TRUE, silent=FALSE, callFrom=NULL) {
  ## was previously cumulCountPlot 
  ##
  argN <- deparse(substitute(roc))
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="staggerdCountsPlot")
  msg <- " argument 'roc' should be matrix (or data.frame) with threshold- and count data (min 2 columns and 2 rows) !"
  if(length(dim(roc)) !=2) stop(msg) else if(any(dim(roc) <2)) stop(msg)
  ## look for counting data
  if(length(countsCol) >0) {
    countsCol <- wrMisc::extrColsDeX(roc,extrCol=countsCol,doExtractCols=FALSE,callFrom=fxNa,silent=silent) 
  } else {
    chColCo <- grep(paste("^",fixedCountPat,sep=""),colnames(roc)) 
    countsCol <- grep(if(length(stats::na.omit(chColCo)) >0) paste("^",fixedCountPat,varCountNa,sep="") else varCountNa, colnames(roc)) }    
  ## look for column with thresholds   threColumn
  if(length(threColumn) <1) stop("argument 'threColumn' seems to be empty")
  if(length(threColumn) >1) threColumn <- threColumn[1]
  threColumn <- wrMisc::extrColsDeX(roc,extrCol=threColumn,doExtractCols=FALSE,callFrom=fxNa,silent=silent) 
  if(length(varCountNa) <1) varCountNa <- colnames(roc)[countsCol]
  ## sort
  newOrd <- order(roc[,threColumn],decreasing=!sortAscending)
  if(!identical(1:nrow(roc),newOrd)) roc <- roc[newOrd,]
  ## eliminiate all 0 count line
  chAll0 <- rowSums(roc[,countsCol]) <= 0
  if(any(chAll0)) { if(!silent & sum(!chAll0) < nrow(roc)/3) message(fxNa, "reducing from ",nrow(roc)," -> ",sum(!chAll0)," rows")
    roc <- roc[which(!chAll0),] }
  maxCoA <- max(rowSums(roc[,countsCol],na.rm=TRUE))
  alph <- roc[,threColumn]
  if(is.null(col)) col <- if(length(varCountNa)==3) {grDevices::rgb(red=c(243,165,240),green=c(184,154,240),blue=c(107,198,150),maxColorValue=255) # pale orange (Ec), purple (Sc),  yellow (Hs)
    } else 1+(1:length(countsCol)) 
  maxCo1 <- max(roc[,countsCol[1]],na.rm=TRUE)
  revInd <- c(1:nrow(roc),nrow(roc):1)
  al2 <- roc[revInd,threColumn]
  if(is.null(tit)) tit <- paste("species counts ",wrMisc::pasteC(varCountNa,quoteC="'")," of ",argN)
  ## help for log-scale
  tmp <- max(0.0003, min(roc[,threColumn]))
  graphics::plot(range(alph),c(0,maxCoA),type="n",xlab=paste("cutoff ",names(threColumn)),ylab="count",main=tit,xaxs="i",yaxs="i", if(logScale) xlim=c(tmp,max(alph,na.rm=TRUE)), log= if(logScale) "x" else "",las=1)
  graphics::polygon( c(alph,range(alph)[c(2,2,1)]), c(roc[,countsCol[1]],maxCo1,0,0),col=col[1],border=grDevices::grey(0.6))        # OK
  if(length(countsCol) >1) for(i in 2:length(countsCol)) graphics::polygon( al2, c(if(i >2) rowSums(roc[,countsCol[(i-1):1]]) else roc[,countsCol[1]],
    rowSums(roc[nrow(roc):1,countsCol[i:1]])),col=col[i])
  if(length(vertLine) >0) {graphics::abline(v=0.05,lty=2,col=grDevices::grey(0.4))           # display alpha
    graphics::mtext(paste(names(threColumn),"=",vertLine),at=vertLine,side=3,line=if(logScale) -1 else 0.1,cex=0.7,las=las.alph)}
  if(displMaxSpec) { z <- roc[nrow(roc),countsCol]
    z <- rbind(z,cumsum(z))
    graphics::mtext(paste("max",varCountNa),side=4,at=z[2,]-z[1,]/2, cex=0.6,las=3,line=0.1)
    graphics::mtext(roc[nrow(roc),countsCol],side=4,at=z[2,]-z[1,]/2, cex=0.7,las=3,line=-1)}
  graphics::legend("topleft",if(length(varCountNa)==length(countsCol)) varCountNa else colnames(roc)[countsCol], text.col=1,pch=22,col=1,pt.bg=col,cex=0.9,pt.cex=1.6,xjust=0.5,yjust=0.5)        
  }
      
