#' Plot linear regression and confidence interval of regression 
#' 
#' This function provides help to display a series of bivariate points given in 'dat' (multiple data formats possible), to model a linear regression and plot the results.
#' Furthermore, a confidence interval to the regression may be added to the plot, regression parameters get be displayed.
#'
#' @param dat (numeric, data.frame or list) main data to plot/inspect. If numeric 'dat' will be used as dependent variable (y-data)
#'  together with numeric 'indepVarLst' (independent variable); if list, then list-elments \code{indepVarLst} and \code{dependVar} will be used; if matrix, the the 1st and 2nd colum will be used
#' @param indepVarLst (character) if 'dat' is list, this designes the list element with the explanatory or independent variable (ie the variable used for explaining, typically x-data) 
#' @param dependVar (character) if 'dat' is list, this designes the list element with dependent variable (ie the variable to be explained, typically y-data) to test
#' @param cusTxt (character) optional custom text to display in subtitle (instead of p-value to H0: slope.regression=0)
#' @param regrLty (integer) line type for regression
#' @param regrLwd (integer) line width for regression
#' @param regrCol (integer) color of regression-line
#' @param confInt (numeric, between 0 and 1) the probabiity alpha for the regression interval, if \code{NULL} no confidence intervall will be plotted/calculated
#' @param confCol (character) (background) color for confidence-interval
#' @param xLab (character) optional custom x-label
#' @param yLab (character) optional custom y-label
#' @param xLim (numeric) custom limit for x-axis (see also \code{\link[graphics]{par}})
#' @param yLim (numeric) custom limit for y-axis (see also \code{\link[graphics]{par}})
#' @param tit (character) optional title
#' @param nSignif (integer) number of significant digits for regression parameters in subtitle of plot
#' @param col (integer or character) custom color for points (choose \code{NULL} for not plotting the actual data)
#' @param pch (integer or character) type of symbol for points (see also \code{\link[graphics]{par}})
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This functions simply plots (to the current graphical devce); an invisible list containing $data, $linRegr, $confInterval (if calculated) may be returned, too
#' @seealso \code{\link[wrMisc]{exclExtrValues}} for decision of potential outliers; \code{\link[graphics]{hist}}, \code{\link{vioplotW}}
#' @examples
#' set.seed(2020); dat1 <- rep(1:6,each=2) +runif(12,0,1)
#' plotLinReg(dat1, gl(6,2))
#' ## extract elements out of list :
#' li2 <- list(aa=gl(5,2), bb=dat1[1:10])
#' plotLinReg(li2, indepVarLst="aa", dependVar="bb")
#' @export  
plotLinReg <- function(dat, indepVarLst=NULL, dependVar=NULL, cusTxt=NULL, regrLty=1, regrLwd=1, regrCol=1, confInt=0.95,
  confCol=NULL, xLab=NULL, yLab=NULL, xLim=NULL, yLim=NULL, tit=NULL, nSignif=3, col=1, pch=1, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## plot linear regression for single gene/protein based on list containing data & annotation for multiple proteins/genes/elements
  ##
  argNa <- c(deparse(substitute(dat)), deparse(substitute(indepVarLst)), deparse(substitute(dependVar)))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotLinReg")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  opar <- graphics::par(no.readonly=TRUE)

  asNumDf <- function(x, colNa=c("x","y")) {   
    ## check matrix or data.frame with 2 columns if numeric, try to convert to data.frame of 2 numeric
    if(!is.data.frame(x)) x <- as.data.frame(x[,1:2], stringsAsFactors=FALSE)
    if(length(colNa) !=2) stop("Argument 'colNa' must be of length=2")
    chNum <- c(is.numeric(x[,1]), is.numeric(x[,2]))     # check for factors
    if(any(!chNum)) for(i in which(!chNum)) {
      num <- try(wrMisc::convToNum(x[,i], spaceRemove=TRUE, remove=NULL, silent=silent, callFrom=fxNa), silent=TRUE)
      if("character" %in% class(num)) {
        num <- as.numeric(as.character(as.factor(x[,i])))
        warning("Trouble converting column no ",i," to numeric (",wrMisc::pasteC(utils::head(x[,i])),",  interpreted as ",wrMisc::pasteC(utils::head(num)),")") } 
      x[,i] <- num }
    colnames(x) <- colNa  
    x }
  extrFromList <- function(x, yy, zz, colNa=c("x","y")) {
    ## look for list-elements names 'yy' & 'zz', use their content ('yy' as 'x' and 'zz' as 'y')
    if(all(is.integer(c(yy[1],zz[1])))) {
      if(any(c(yy[1],zz[1]) <1) || length(x) < max(c(yy[1],zz[1]))) stop(" index values for list-elements of 'x' out of range")
    } else {
      ## thus yy& zz are considered/tested as names of x
      msg <- "both arguments 'yy' and 'zz' must correspond to list-elements of 'x'"
      if(length(yy) <1) stop(msg)
      if(length(zz) <1) {    # if 'zz' NULL, try to extract 1st & 2nd col of x$yy
        isBad <- TRUE
        if(length(dim(x[[yy]])) >1) if(ncol(x[[yy]]) >1) { isBad <- FALSE
          x <- as.data.frame(x[[yy]][,1:2], stringsAsFactors=FALSE) } 
        if(isBad) stop(msg)
      } else {      
        chNa1 <- c(yy[1], zz[1]) %in% names(x)
        if(!all(chNa1)) stop("Cannot find ",wrMisc::pasteC(c(yy[1],zz[1])[which(!chNa1)], quoteC="'")," in list 'x'")
        x <- data.frame(x[[yy]], x[[zz]], stringsAsFactors=FALSE)}
    }    
    colnames(x) <- colNa
    x }    
  extrFromMatr <- function(x, yy, zz, name1=c("y","ordinate","dat","measure","pred","depend"), name2=c("x","abscissa","grp","grp2","dat2","obs","indep"),
    colNa=c("x","y"), silent=silent, fxNa=fxNa) {
    ## extract column specified in 'zz' or look among names in 'name1'
    ## if no name found for 'yy' just avoid using using the column found for 'zz'
    chColNa1 <- wrMisc::extrColsDeX(x, extrCol=list(if(is.null(zz)) name1 else zz), doExtractCols=FALSE, callFrom=fxNa,silent=silent)
    chColNa2 <- wrMisc::extrColsDeX(x, extrCol=list(if(is.null(yy)) name2 else yy), doExtractCols=FALSE, callFrom=fxNa,silent=silent)
    if(!all(chColNa1,chColNa2)) stop("Cannot find column-names to use from 'x'")
    x <- data.frame(x[,if(length(chColNa2) >0) chColNa2[1] else { if(chColNa1[1]==2) 1 else chColNa1[1]+1}], x[,chColNa1[2]],stringsAsFactors=FALSE) 
    colnames(x) <- colNa
    x }
  ##      
  msg <- df0<- NULL                               # initialize
  ## check main input : see if matrix providing x & y
  if(length(dat) <1) { 
    msg <- c(" incomplete data, nothing to do")
  } else { 
    if(!is.list(dat) && length(dat) >2 && length(indepVarLst) >2 && length(dim(dat)) <1 && length(dim(indepVarLst)) <1) {
      ## simplest case: both x & y as separate numeric vector (of same length)
      if(length(dat) !=length(indepVarLst)) stop("Length of 'dat' and 'indepVarLst' don't match !")
      df0 <- data.frame(x=indepVarLst, y=dat, stringsAsFactors=FALSE) 
      argNa[4:5] <- argNa[2:1]
    } else {
      ## if S3 (from limma) or (other) list, need to locate-list-elements
      if(is.list(dat)) {
        df0 <- extrFromList(dat, indepVarLst, dependVar)
        argNa[4:5] <- argNa[2:3]
      } else {
        ## dat is not a list, assume dat is matrix
        if(length(dim(dat)) >1) { if(ncol(dat) >1) {
          df0 <- extrFromMatr(dat, indepVarLst, dependVar, silent,fxNa)
        } else {
          ## dat is neiter a matrix, try extracting 1st col of indepVarLst
          df0 <- data.frame(x=if(length(dim(indepVarLst)) >1) indepVarLst[,1] else indepVarLst, y=dat, stringsAsFactors=FALSE) }
        argNa[4:5] <- argNa[2:3]
        } else msg <- "unknown format of 'dat'" }}}
  ## check if plot can be produced
  if(length(msg) >0 || length(df0) <1) message(fxNa,"Can't plot",msg) else {
  ## make linear model
  df0 <- asNumDf(df0)
  lm0 <- stats::lm(y ~ x, data=df0) 
  if(length(lm0$coefficients) >2) message(fxNa,"Bizzare : The regression model was expected as 2 coefficients, but has ",
    length(lm0$coefficients)," coefficients ",wrMisc::pasteC(names(lm0$coefficients),quoteC="'"))
  ## start plotting
  argNa[4:5] <- gsub("\"","",argNa[4:5])                   # remove protected 'double' quotes
  if(is.null(xLab)) xLab <- if(argNa[4]=="NULL") "x" else argNa[4]                 # explanatory variable
  if(is.null(yLab)) yLab <- if(argNa[5]=="NULL" | argNa[5]==xLab) "y" else argNa[5]  
  tmp <- try(graphics::plot(y ~ x, data=df0, las=1, xlab=xLab, ylab=yLab, pch=pch,col=col, main=tit), silent=TRUE)
  if(inherits(tmp, "try-error"))  warning(fxNa," Plot cannot be produced") else {
    graphics::abline(lm0, lty=regrLty, lwd=regrLwd, col=regrCol)
    suplTx <- paste(c("; ",if(length(cusTxt) <1) paste("p.slope =",signif(stats::coef(summary(lm0))[2,"Pr(>|t|)"],2)) else cusTxt), collapse=" ")
    graphics::mtext(paste("regression (rounded): y =",signif(stats::coef(lm0)[2],nSignif)," x +",signif(stats::coef(lm0)[1],nSignif),suplTx,
      ",  r2=",signif(stats::cor(df0$y,df0$x)^2,nSignif)),cex=0.75,line=0.15)
    if(length(confInt) >0) { ra <- c(range(df0$x,na.rm=TRUE), abs(mean(df0$x,na.rm=TRUE)))
      newx <- seq(ra[1]-0.05*ra[3],ra[2]+0.05*ra[3],length.out=200)
      if(length(confCol) <1) confCol <- grDevices::rgb(0.3,0.3,0.3,0.07)           # (background) color for confidence-interval
      confInterval <- stats::predict(lm0, newdata=data.frame(x=newx), interval="confidence", level=confInt)   # can do single conf interv at a time ..
      graphics::polygon(cbind(x=c(newx,rev(newx)),y=c(confInterval[,"lwr"],confInterval[length(newx):1,"upr"])),col=confCol,border=NA)
      graphics::points(y ~ x, df0, col=col)
      graphics::mtext(paste("  confidence interval at ",100*confInt,"% shown"), line=-1.05,cex=0.65,adj=0,col=wrMisc::convColorToTransp(confCol,alph=240))  } }
  }
  invisible(list(data=df0,linRegr=lm0,if(length(confInt) >0) confInterval=confInterval)) }
  
