#' Display numeric content of matrix as image  
#'
#' To get a quick overview of the distribution of data and, in particular, of local phenomena it is useful to express numeric values as colored boxes.
#' Such an output can be referred to as heatmap.
#' The function \code{\link[graphics]{image}} provides the basic support to do so.  
#' To do this more conveniently, the function \code{imageW} offers additional options for displaying row- and column-names or displaying NA-values as custom-color.
#'
#' @param dat (list) main input
#' @param col (character or integer) colors, default is heat.colors in 15 shades
#' @param rowNa (character) optional custom rownames
#' @param colNa (character) optional custom colnames
#' @param tit (character) custom figure title
#' @param cexXlab (numeric) cex-like expansion factor for x-axis labels  (see also \code{\link[graphics]{par}})
#' @param cexYlab (numeric) cex-like expansion factor for y-axis labels  (see also \code{\link[graphics]{par}})
#' @param cexTit (numeric) cex-like expansion factor for title  (see also \code{\link[graphics]{par}})
#' @param NAcol (character or integer) custom color fro NA-values, default is grey
#' @param las (numeric) style of axis labels (see also \code{\link[graphics]{par}})
#' @seealso \code{\link[graphics]{image}}, heatmaps including hierarchical clustering \code{\link[stats]{heatmap}} or \code{heatmap.2} from  \href{https://CRAN.R-project.org/package=gplots}{gplots}   
#' @return graphical output only
#' @examples
#' imageW(as.matrix(iris[1:40,1:4]))
#' @export
imageW <- function(dat, col=NULL, rowNa=NULL, colNa=NULL, tit=NULL, cexXlab=0.7,cexYlab=1, cexTit=1.6, NAcol=grDevices::grey(0.8),las=2) {
  ## improved version if image()
  dat <- dat[,ncol(dat):1]
  if(length(col) <1) col <- rev(grDevices::heat.colors(15))
  if(any(is.na(dat)) & length(NAcol) >0) { if(!is.matrix(dat)) dat <- as.matrix(dat)
    mi <- min(dat,na.rm=TRUE)
    dat[which(is.na(dat))] <- min(dat,na.rm=TRUE) -diff(range(dat,na.rm=TRUE))*1.1/(length(col))
     col <- c(rep(NAcol,1),col) }
  if(length(rowNa) <nrow(dat)) rowNa <- rownames(dat)
  if(length(rowNa) <1) rowNa <- 1:nrow(dat)
  if(length(colNa) <ncol(dat)) colNa <- colnames(dat)
  if(length(colNa) <1) rowNa <- 1:ncol(dat)  
  graphics::image(dat, col=col, xaxt="n", yaxt="n", main=tit, cex.main=cexTit)
  graphics::mtext(at=(0:(length(colNa)-1))/(length(colNa)-1), colNa, cex=cexYlab, side=2, line=0.3, las=las)   # on left
  graphics::mtext(at=(0:(length(rowNa)-1))/(length(rowNa)-1), rowNa, cex=cexXlab, side=1, line=0.3, las=las)   # on bottom
  graphics::box(col=grDevices::grey(0.8)) }
  
