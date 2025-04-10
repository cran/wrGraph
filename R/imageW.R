#' Display Numeric Content Of Matrix As Image  
#'
#' To get a quick overview of the spatial distribution of smaller data-sets it may be useful to display numeric values as colored boxes.
#' Such an output may also be referred to as heatmap (note that the term 'heatmap' is frequently associated with graphical display of hierarchcal clustering results).
#' The function \code{\link[graphics]{image}} provides the basic support to do so (ie heatmap without rearranging rows and columns by clustering).  
#' To do this more conveniently, the function \code{imageW} offers additional options for displaying row- and column-names or displaying NA-values as custom-color.
#'
#' @details 
#' This function allows two modes of operation : 1) plotting using standard R -graphics or 2) using the framework of grid- and lattice-graphics (since version 1.2.6). 
#' The latter version allows integrating a legend for the color-scale and adding grid-lines, rotation of axis-labels or removing tick-marks.
#' Please note that sometimes the center-color segment may not end up directly with the center of the scale, in this case you may adjust using the argument \code{centColShift=-1}
#'
#' @param data (matrix or data.frame) main input
#' @param latticeVersion (logical) use lattice for plotting (this will include a color-legend)
#' @param transp (logical) decide if data should get transposed (if \code{TRUE} the data will be displayed exacetly same order as when printing the values as table); 
#'  set to \code{FALSE} to get behaviour prior to version 1.3.0. 
#' @param col (character or integer) colors; in lattice version 2 or 3 color-names to define central- and end-points of gradient (starting with color for lowest values, optional central color and color for highest values), default is 60 shades 'RdYlBu' RColorBrewer, if 'heat.colors' use  heat.colors in min 15 shades
#' @param NAcol (character or integer) custom color for NA-values, default is light grey
#' @param rowNa (character) optional custom rownames
#' @param colNa (character) optional custom colnames
#' @param tit (character) custom figure title
#' @param xLabVal (character) optional custom text for x-axis 'values' (multiple values/names instead of counters, replaces argument \code{xLab} in older versions)
#' @param yLabVal (character) optional custom text for y-axis 'values' (multiple values/names instead of counters, replaces argument \code{yLab} in older versions)
#' @param xLab (character, length=1) optional custom text for x-axis label (so far fixed color & fontsiez)
#' @param yLab (character, length=1) optional custom text for y-axis label
#' @param las (numeric) style of axis labels (see also \code{\link[graphics]{par}}); in case of \code{latticeVersion=TRUE} this argument will override default \code{rotXlab=0} and/or  \code{rotYlab=0}
#' @param nColor (integer, only used in lattice version) number of color-blocks in color gradient (made based on central- and end-points from \code{col} 
#' @param balanceCol (logical, only used in lattice version) if \code{TRUE} the color-radient aims to color the value closest to 0 with the center color (from \code{col} (default gray)
#' @param gridCol (character, only used in lattice version) define color of grid 
#' @param gridLty (integer, only used in lattice version) define line-type of grid (see also lty \code{\link[graphics]{par}})
#' @param centColShift (integer, only used in lattice version) shift central (default grey) color element for negative scale up or down (ie increase or reduce number of color-blocks for negatve values), 
#'   used for correcting automatic scaling rounding issues to ensure the central elements captures 0  
#' @param cexDispl (numeric, length=1, only used in lattice version) define cex size for displaying (rounded) values in plot, set to \code{NULL} for omitting
#' @param panel.background.col (character, only used in lattice version)
#' @param supLat (list, only used in lattice version) additional arguments/parameters passed to \code{levelplot} - currently not activated  
#' @param rotXlab (numeric, 0 - 360, lattice version only) control rotation of x-axis values
#' @param rotYlab (numeric, 0 - 360, lattice version only) control rotation of y-axis values
#' @param cexTit (numeric) cex-like expansion factor for title  (see also \code{\link[graphics]{par}})
#' @param cexAxs (numeric) cex-like expansion factor for x- and y-axis text/labels (see also \code{\link[graphics]{par}})
#' @param cexXlab (numeric) cex-like expansion factor for x-axis labels  (see also \code{\link[graphics]{par}})
#' @param cexYlab (numeric) cex-like expansion factor for y-axis labels  (see also \code{\link[graphics]{par}})
#' @param showValues (logical or numeric) optional display of values from data, if contains eg \code{cex=0.5} this value will be used as expansion factor 
#'  (otherwise default to 0.6); if contains eg \code{digits=3} this value displayed will be rounded to this number of significant digits  (otherwise default to 4)
#' 
#' @param Xtck (numeric or logical) expansion factor for length of tick-marks on x-axis (default=0 for no tick-marks)
#' @param Ytck (numeric or logical) expansion factor for length of tick-marks on y-axis
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' 
#' @seealso \code{\link[graphics]{image}}, for the lattice version \code{\link[lattice]{levelplot}}, heatmaps including hierarchical clustering \code{\link[stats]{heatmap}} or \code{heatmap.2} from package \href{https://CRAN.R-project.org/package=gplots}{gplots}   
#' @return This function plots in image (to the current graphical device) as \code{image} does
#' @examples
#' imageW(iris[1:40,1:4], transp=FALSE, tit="Iris (head)")
#' imageW(iris[1:20,1:4], latticeVersion=TRUE, col=c("blue","red"), 
#'   rotXlab=45, yLab="Observation no", tit="Iris (head)")
#' @export
imageW <- function(data, latticeVersion=FALSE, transp=TRUE, NAcol="grey95", tit=NULL, rowNa=NULL, colNa=NULL, xLab=NULL, yLab=NULL, xLabVal=NULL, yLabVal=NULL, las=2, 
  col=NULL, nColor=9, balanceCol=TRUE, gridCol="grey75", gridLty=1, centColShift=0, cexDispl=NULL, panel.background.col="white", supLat=list(),
  rotXlab=0, rotYlab=0, cexTit=1.6, cexAxs=NULL, cexXlab=0.7, cexYlab=0.9, showValues=FALSE, Xtck=0, Ytck=0, silent=FALSE, debug=FALSE, callFrom=NULL) { 
  ## improved version if image() or  levelplot()
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="imageW")
  argNa <- deparse(substitute(data))
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  doPlot <- length(data) > 0 && any(c("numeric","matrix","data.frame") %in% class(data))
  if(doPlot) {if(length(dim(data)) <2) data <- as.matrix(data, ncol=1, dimnames=list(names(data), NULL))
    chNum <- is.numeric(data[1,1])
    if(!chNum) {iniCl <- class(data); iniDimNa <- dimnames(data)
      data <- try(as.numeric(as.matrix), silent=TRUE)
      if(inherits(data, "try-error")) doPlot <- FALSE else {
        ##  managed to convert to numeric successfully
        data <- matrix(data, nrow=length(iniDimNa[[1]]), dimnames=iniDimNa)
        if("data.frame" %in% iniCl) data <- as.data.frame(data) }
    }
  }  

  transp <- !isFALSE(transp)
#  if(length(dim(data)) <2) data <- try(matrix(as.numeric(data), ncol=1, dimnames=list(names(data), NULL)), silent=TRUE)
#  if(inherits(data, "try-error")) doPlot <- FALSE else {
#    if(is.data.frame(data) && doPlot) {doPlot <- is.numeric(as.matrix(data)); data <- as.matrix(data)}}
  if(debug) {message(fxNa," xiW0   debug: ",debug); xiW0 <- list(data=data,xLabVal=xLabVal,yLabVal=yLabVal,transp=transp,latticeVersion=latticeVersion,doPlot=doPlot, tit=tit,rowNa=rowNa,colNa=colNa,showValues=showValues)}
  if(doPlot) {    
    ## checks & adjust
    if(length(xLabVal) >1 && !all(is.na(xLabVal))) {
      if(transp) { 
        if(length(xLabVal) ==ncol(data) && length(colNa) <1 ) { colNa <- xLabVal; xLabVal <- NA
          if(!silent) message(fxNa,"(88) It seems you meant 'colNa' when using argument 'xLabVal' (interpreting as such) ...") 
        } else { if(length(xLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'xLabVal'"); xLabVal <- NA }
      } else {   
        if(length(xLabVal) ==nrow(data) && length(rowNa) <1 ) { rowNa <- xLabVal; xLabVal <- NA
          if(!silent) message(fxNa,"(92) It seems you meant 'rowNa' when using argument 'xLabVal' (interpreting as such) ...") 
        } else { if(length(xLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'xLabVal'"); xLabVal <- NA }
      } 
    } else { if(length(xLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'xLabVal'"); xLabVal <- NA }
    
    if(length(yLabVal) >1 && !all(is.na(yLabVal))) {
      if(transp) { 
        if(length(yLabVal) ==nrow(data) && length(rowNa) <1 ) { rowNa <- yLabVal; yLabVal <- NA
          if(!silent) message(fxNa,"(100) It seems you meant 'rowNa' when using argument 'yLabVal' (interpreting as such) ...")
        } else { if(length(yLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'yLabVal'"); yLabVal <- NA }
      } else {  
        if(length(yLabVal) ==ncol(data) && length(colNa) <1 ) { colNa <- yLabVal; yLabVal <- NA
          if(!silent) message(fxNa,"(104) It seems you meant 'colNa' when using argument 'yLabVal' (interpreting as such) ...")
        } else { if(length(yLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'xLabVal'"); yLabVal <- NA}           
      }      
    } else { if(length(yLabVal) >1) if(!silent) message(fxNa,"Invalid entry for 'xLabVal'"); yLabVal <- NA}
                                                                                                                     
    if(length(rowNa) < nrow(data)) rowNa <- if(nrow(data) >1) {if(length(rownames(data)) >0) rownames(data) else 1:nrow(data)} else ""
    if(length(colNa) < ncol(data)) colNa <- if(ncol(data) >1) {if(length(colnames(data)) >0) colnames(data) else 1:ncol(data)} else ""
    if(debug) {message(fxNa," xiW1"); xiW1 <- list(data=data,xLabVal=xLabVal,yLabVal=yLabVal,transp=transp,latticeVersion=latticeVersion, tit=tit,rowNa=rowNa,colNa=colNa,showValues=showValues)}
    if(isTRUE(latticeVersion)) {
      ## reformat input
      if(!transp) data <- t(data)                 #  was initially written for transp=T, re-transform if not chosen
      if(length(rotXlab)==0 && any(las %in% c(2,3))) rotXlab <- 0
      if(length(rotYlab)==0 && any(las %in% c(0,3))) rotYlab <- 0
      ma2 <- expand.grid(1:ncol(data), 1:nrow(data))
      ma2 <- cbind(ma2, as.numeric(t(data[nrow(data):1,])))
      colnames(ma2) <- c("x","y","z")
      if(any(is.na(xLabVal))) xLabVal <- NULL
      if(any(is.na(yLabVal))) yLabVal <- NULL
      ## colors
      if(length(col) <2) col <- c("blue","grey80","red")
      nCol2 <- try(round(nColor[1]), silent=TRUE)
      msg <- "Note: Argument 'nColor' should contain integer at least as high as numbers of colors defined to pass through; resetting to default=9"
      if(inherits(nCol2, "try-error")) { message(fxNa,msg); nCol2 <- 9 }
      if(nCol2 < length(col)) { message(fxNa,msg); nCol2 <- 9 }
      if(debug) {message(fxNa," xiW2"); xiW2 <- list(data=data,xLabVal=xLabVal,yLabVal=yLabVal,transp=transp,latticeVersion=latticeVersion, nCol2=nCol2,tit=tit,rowNa=rowNa,colNa=colNa,showValues=showValues)}

      miMa <- range(as.matrix(data), na.rm=TRUE)
      width <- (miMa[2] - miMa[1])/ nCol2
      bre <- miMa[1] + (0:nCol2) *width           # breaks
      clo0 <- which.min(abs(as.numeric(as.matrix(data))))    # (first) value closest to 0, try to include in grey segm
      clo0br <- min(which(bre >= as.numeric(as.matrix(data))[clo0]))   #+ (-1:0)  # upper break/bound for center color (close/including 0)
      if(debug) {message(fxNa," xiW2b")}

      if(clo0br >1 && clo0br < length(bre) -1 && length(col) >2) {   # some values in lower & upper gradient
        maxLe <- max(clo0br -1, length(bre) -clo0br) -as.integer(balanceCol)  
        negCol <- try(grDevices::colorRampPalette(col[1:2])(if(balanceCol) maxLe else length(clo0br -1)), silent=TRUE)
        if(inherits(negCol, "try-error")) { negCol <- NULL
          if(!silent) message(fxNa,"Invalid color-gradient for neg values")
        }
        negCol <- negCol[-length(negCol)]                              # max neg-col -> grey (wo defined grey);  remove 'grey' from last position
        posCol <- try((grDevices::colorRampPalette(col[2:3])(if(balanceCol) maxLe else length(length(bre) -1 -clo0br))), silent=TRUE)  # (grey -> max pos-col)
        if(inherits(posCol, "try-error")) { 
          if(!silent) warning(fxNa,"Invalid color-gradient for pos values")
        }
        if(debug) {message(fxNa, "/1 clo0br ",clo0br,"   max nCol ",nCol2,"   le negCol ",length(negCol),"   le posCol ",length(posCol)," xiW2c")}
        if(balanceCol) {
          centColShift <- if(length(centColShift) <1 || !is.numeric(centColShift)) 0 else as.integer(centColShift)
          .keepLastN <- function(x,lastN) x[(length(x) -lastN +1):length(x)]
          if(length(negCol) != clo0br -2 +centColShift) {
            if(debug) message(fxNa,"Correct negCol (prev=",length(negCol),") centColShift=",centColShift," to : ",clo0br -2 +centColShift)
            if(length(negCol) > clo0br -2 +centColShift) negCol <- .keepLastN(negCol, clo0br -2 +centColShift)
            if(length(negCol) < clo0br -2 +centColShift) {negCol <- grDevices::colorRampPalette(col[1:2])(clo0br -1 +centColShift)
              negCol <- negCol[-length(negCol)] }
          }
          if(length(posCol) != length(bre) -length(negCol) -1) {
            if(debug) message(fxNa,"Corr posCol (prev ",length(posCol),") to ",maxLe + centColShift," to ",length(bre) -length(negCol) -1)
            if(length(posCol) > length(bre) -length(negCol) -1) posCol <- posCol[1:(length(bre) -clo0br)]
            if(length(posCol) < length(bre) -length(negCol) -1) {
              posCol <- grDevices::colorRampPalette(col[2:3])(length(bre) -length(negCol) -1) }
          }
        } 
        cols <- c(negCol, posCol)
        if(debug) message(fxNa, "/2 clo0br ",clo0br,"   max nCol ",nCol2,"  le cols ",length(cols),"   le negCol ",length(negCol),"   le posCol ",length(posCol))
      } else {   # plain color gradient
         cols <- if(length(col)==2) grDevices::colorRampPalette(col[1:2])(length(bre) -1) else {
           c(grDevices::colorRampPalette(col[1:2])(floor(length(bre)/2)), (grDevices::colorRampPalette(col[2:3])(length(bre) -floor(length(bre)/2)))[-1])
        }
      }
      ##  
      myPanel <- function(...) {
        grid::grid.rect(gp=grid::gpar(col=NA, fill=NAcol))       # fill NA
        lattice::panel.levelplot(...)
        argXYZ <- list(...)
        if(length(cexDispl)==1 && is.numeric(cexDispl)) lattice::panel.text(argXYZ$x, argXYZ$y, signif(argXYZ$z,2), cex=cexDispl)      # add rounded numeric value
        if(any(is.na(gridCol))) gridCol <- NULL
        chGri <- (1:6) %in% gridLty
        if(length(gridCol) >0 && any(chGri)) {                # add grid-lines
          lattice::panel.abline(h=0.5 +1:(nrow(data) -1), col=gridCol, lty=gridLty)     # vertical 
          lattice::panel.abline(v=0.5 +1:(ncol(data) -1), col=gridCol, lty=gridLty) }   # hor
        ## add axis label (new oct24)
        if(length(xLab)==1 && !is.na(xLab)) lattice::panel.text(x=mean(argXYZ$x), y=0.02, xLab, gp=grid::gpar(fontsize=14), just="center", check.overlap=TRUE) 
        if(length(yLab)==1 && !is.na(yLab)) lattice::panel.text(x=0.02, mean(argXYZ$y), yLab, gp=grid::gpar(fontsize=14), vjust=0.5, rot=90, check.overlap=TRUE)  # text may be too much to left
      } 
          #panel.text(x=mytable$Xcoord, y=mytable$Ycoord, mytable$Labels)              
      if(debug) {message(fxNa," xiW3"); xiW3 <- list(data=data,xLabVal=xLabVal,yLabVal=yLabVal,transp=transp,latticeVersion=latticeVersion, tit=tit,rowNa=rowNa,colNa=colNa,showValues=showValues,ma2=ma2,cols=cols,cexXlab=cexXlab,rotXlab=rotXlab,Xtck=Xtck,Ytck=Ytck)}

      ## lattice levelplot
      plo <- try(lattice::levelplot(z ~ x *y, data=ma2, aspect=nrow(data)/ncol(data), col.regions=cols,
        region=TRUE, cuts=length(cols) -1, xlab=yLabVal, ylab=xLabVal, main=tit,
        scales=list(relation="free", x=list(at=1:ncol(data), labels=if(transp) colNa else rowNa, cex=cexXlab, rot=rotXlab, tck=as.numeric(Xtck)), 
          y=list(at=nrow(data):1, labels=if(transp) rowNa else colNa, cex=cexYlab, rot=rotYlab, tck=as.numeric(Ytck))),  # axis labels 
        par.settings=list(axis.line=list(col='black')), 
        panel=myPanel ), silent=TRUE)
      if(inherits(plo, "try-error")) message(fxNa, "Failed to plot !!\n", plo) else {                 # new oct24
        #if(length(xLab)==1 && !is.na(xLab)) plo <- plo + grid::grid.text(xLab, y=0.02, gp=grid::gpar(fontsize=14), just="center", check.overlap=TRUE) 
        #if(length(yLab)==1 && !is.na(yLab)) plo <- plo + grid::grid.text(yLab, x=0.02, gp=grid::gpar(fontsize=14), vjust=0.5, rot=90, check.overlap=TRUE)  # text may be too much to left
      }  
         
    } else {
      ## (until v1.2.5) standard graphics version  (ie non-lattice)
      if(debug) message(fxNa," (175)  xVal=",xLabVal,"   yVal=",yLabVal)
      if(transp) data <- t(data)
      if(ncol(data) >1) data <- data[,ncol(data):1]                         # reverse for intuitive display left -> right
      if(identical(col,"heat.colors") || identical(col,"heatColors")) col <- rev(grDevices::heat.colors(sort(c(15, prod(dim(data)) +2))[2] ))
      chRCo <- requireNamespace("RColorBrewer", quietly=TRUE) 
      msgRCo <- c(fxNa,"Package 'RColorBrewer' not installed",", ignore argument 'col'")
      if(identical(col,"YlOrRd"))  {if(chRCo) col <- RColorBrewer::brewer.pal(9,"YlOrRd") else { col <- NULL; if(!silent) message(msgRCo) }}
      if(identical(col,"RdYlGn"))  {if(chRCo) col <- RColorBrewer::brewer.pal(11,"RdYlGn") else { col <- NULL; if(!silent) message(msgRCo) }}
      if(identical(col,"Spectral"))  {if(chRCo) col <- RColorBrewer::brewer.pal(11,"Spectral") else { col <- NULL; if(!silent) message(msgRCo) }}
      if(identical(col,"RdBu"))  {if(chRCo) col <- RColorBrewer::brewer.pal(11,"RdBu") else { col <- NULL; if(!silent) message(msgRCo) }}  
      if(length(col) <1) { if(!chRCo) message(msgRCo[1:2],"Using rainbow colors instead of 'RdYlBu'") 
        col <- if(chRCo) grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(60) else grDevices::rainbow(60)}  
      chNa <- is.na(data)  
      if(any(chNa) && length(NAcol) >0) { if(!is.matrix(data)) data <- as.matrix(data)
        mi <- min(data, na.rm=TRUE)
        ## mark NAs
        if(any(chNa)) data[which(chNa)] <- min(data, na.rm=TRUE) -diff(range(data, na.rm=TRUE))*1.1/(length(col))
        col <- c(NAcol, col) }
      ## main plot
      yAt <- (0:(length(rowNa)-1))/(length(rowNa)-1)
      if(is.data.frame(data)) data <- as.matrix(data)
      if(debug) {message(fxNa," (197)  xVal=",xLabVal,"   yVal=",yLabVal, "  xiW4"); xiW4 <- list(data=data, col=col, xLabVal=xLabVal,yLabVal=yLabVal,transp=transp,tit=tit,rowNa=rowNa,colNa=colNa,yAt=yAt,showValues=showValues,cexTit=cexTit)}
      plo <- try(graphics::image(data, col=col, xaxt="n", yaxt="n", main=tit, cex.main=cexTit, xlab=if(transp) yLabVal else xLabVal, ylab=if(transp) xLabVal else yLabVal), silent=TRUE)
      if(inherits(plo, "try-error")) message(fxNa,"Failed to plot !!\n", plo) else {
        graphics::mtext(at=(0:(length(colNa)-1))/(length(colNa)-1), colNa, side=if(transp) 1 else 2, line=0.3, las=las, cex=cexYlab)   # 'colNames'
        graphics::mtext(at=if(transp) rev(yAt) else yAt, rowNa, side=if(transp) 2 else 1, line=0.3, las=las, cex=cexXlab)           # 'rowNames'
        graphics::box(col=grDevices::grey(0.8)) 
        showVal <- if(length(showValues) >0) !isFALSE(showValues[1]) else FALSE
        if(showVal) {
          nDigNa <- c("nDig","dig","digits","sig","sign","signif")
          nSign <- if(any(nDigNa %in% names(showValues))) showValues[which(names(showValues) %in% nDigNa)[1]] else 4
          cexTxt <- if(any(c("cex") %in% names(showValues))) showValues[which(names(showValues) %in% c("cex"))[1]] else 0.6
          col <- if(any(c("col","color") %in% names(showValues))) showValues[which(names(showValues) %in% c("col","color"))[1]] else 1
          plo <- try(graphics::text( rep(seq(1,0, length.out=ncol(data)), each=nrow(data)), rep(seq(1,0, length.out=nrow(data)), ncol(data)), 
            labels=signif(as.numeric(data), nSign), col=col, cex=cexTxt), silent=TRUE)
          if(inherits(plo, "try-error")) message(fxNa,"Failed adding values to plot !\n", plo) }                
      }       
    }                 # finish base graphics plot
  } else if(!silent) message(fxNa,"Argument 'data' invalid, please furnish matrix or data.frame with min 2 lines & min 1 col") 
}  
               


