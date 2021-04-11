#' Display numeric content of matrix as image  
#'
#' To get a quick overview of the distribution of data and, in particular, of local phenomena it is useful to express numeric values as colored boxes.
#' Such an output can also be referred to as heatmap (note that the term 'hatmap' is also frequently associated with graphical display of hierarchcal clustering results).
#' The function \code{\link[graphics]{image}} provides the basic support to do so (ie heatmap without rearranging rows and columns by clustering).  
#' To do this more conveniently, the function \code{imageW} offers additional options for displaying row- and column-names or displaying NA-values as custom-color.
#'
#' @details 
#' If the main input \code{dat} is numeric vector (an not matrix or data.frame) the values will be displayed as multiple columns and single row.
#'
#' @param dat (matrix or dtaa.frame) main input
#' @param col (character or integer) colors, default is 60 shades 'RdYlBu' RColorBrewer, if 'heat.colors' use  heat.colors in min 15 shades
#' @param rowNa (character) optional custom rownames
#' @param colNa (character) optional custom colnames
#' @param tit (character) custom figure title
#' @param xLab (character) optional custom names for x-axis
#' @param yLab (character) optional custom names for y-axis
#' @param cexXlab (numeric) cex-like expansion factor for x-axis labels  (see also \code{\link[graphics]{par}})
#' @param cexAxs (numeric) cex-like expansion factor for x- and y-axis text/labels (see also \code{\link[graphics]{par}})
#' @param cexYlab (numeric) cex-like expansion factor for y-axis labels  (see also \code{\link[graphics]{par}})
#' @param cexTit (numeric) cex-like expansion factor for title  (see also \code{\link[graphics]{par}})
#' @param NAcol (character or integer) custom color fro NA-values, default is grey
#' @param las (numeric) style of axis labels (see also \code{\link[graphics]{par}})
#' @seealso \code{\link[graphics]{image}}, heatmaps including hierarchical clustering \code{\link[stats]{heatmap}} or \code{heatmap.2} from  \href{https://CRAN.R-project.org/package=gplots}{gplots}   
#' @return graphical output only
#' @examples
#' imageW(as.matrix(iris[1:40,1:4]))
#' @export
imageW <- function(dat, col=NULL, rowNa=NULL, colNa=NULL, tit=NULL, xLab=NA, yLab=NA, cexXlab=0.7, cexAxs=NULL, cexYlab=0.9, cexTit=1.6, NAcol=grDevices::grey(0.8), las=2) {
  ## improved version if image()
  if(length(dim(dat)) <2) dat <- matrix(as.numeric(dat), ncol=1, dimnames=list(names(dat), NULL))
  if(ncol(dat) >1) dat <- dat[,ncol(dat):1]
  if(identical(col,"heat.colors") | identical(col,"heatColors")) col <- rev(grDevices::heat.colors(sort(c(15, prod(dim(dat)) +2))[2] ))
  if(length(col) <1) col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(60)  
  chNa <- is.na(dat)  
  if(any(chNa) & length(NAcol) >0) { if(!is.matrix(dat)) dat <- as.matrix(dat)
    mi <- min(dat,na.rm=TRUE)
    ## mark NAs
    if(any(chNa)) dat[which(chNa)] <- min(dat,na.rm=TRUE) -diff(range(dat,na.rm=TRUE))*1.1/(length(col))
    col <- c(NAcol,col) }
  if(length(dim(dat)) <2) dat <- matrix(dat, ncol=1, dimnames=list(names(dat),NULL))
  if(length(rowNa) <nrow(dat)) rowNa <- rownames(dat) 
  if(length(rowNa) <1) rowNa <- if(length(nrow(dat)) >1) 1:nrow(dat) else ""
  if(length(colNa) <ncol(dat)) colNa <- colnames(dat)
  if(length(colNa) <1) colNa <- if(length(ncol(dat)) >1) 1:ncol(dat) else ""
  if(is.null(xLab)) xLab <- ""
  if(is.null(yLab)) yLab <- ""
  ## main plot
  graphics::image(dat, col=col, xaxt="n", yaxt="n", main=tit, xlab=xLab, ylab=yLab, cex.main=cexTit)
  graphics::mtext(at=(0:(length(colNa)-1))/(length(colNa)-1), colNa, side=2, line=0.3, las=las, cex=cexYlab)   # on left  , cex=cexAxs
  graphics::mtext(at=(0:(length(rowNa)-1))/(length(rowNa)-1), rowNa, side=1, line=0.3, las=las, cex=cexXlab)   # on bottom  , cex=cexAxs
  graphics::box(col=grDevices::grey(0.8)) }
  
