#' MA-plot (differential intensity versus average intensity)   
#'
#' This type of plot is very common in high-throughput biology, see \href{https://en.wikipedia.org/wiki/MA_plot}{MA-plot}.
#' Basically one would like to compare numerous independent measures (ie gene transcript or protein abundance values) of 2 samples/data-sets,
#' it is usual to compare a change ('Minus'=M) versus absolute mean value ('Average'=A).
#' In high-throughput biology data are typically already transformed to log2 and thus, the 'M'-value represents a relative change.
#' Besides, output from statistical testing by \code{\link[wrMisc]{moderTest2grp}} can be directly read to produce MA plots for diagnostic reasons.
#' Please note, that plotting a very number of points in transparency (eg >10000) may take several seconds.
#'	 
#' @param Mvalue (numeric or matrix) data to plot; M-values are typically calculated as difference of log2-abundance values and 'Avalue' the mean of log2-abundance values;
#'   M-values and A-values may be given as 2 columsn of a matrix, in this case the argument \code{Avalue} should remain NULL 
#' @param Avalue (numeric, list or data.frame) if \code{NULL} it is assumed that 2nd column of 'Mvalue' contains the A-values to be used
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied 
#' @param ProjNa (character) custom title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC)
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points as increased shades of gray 
#' @param compNa (character) names of groups compared
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace  
#' @param cexMa (numeric) expansion factor for font-size of title (see also \code{\link[graphics]{par}})
#' @param cexLa (numeric) expansion factor \code{cex} for labels (see also \code{\link[graphics]{par}}) 
#' @param limM (numeric, length=2) range of axis M-values 
#' @param limA (numeric, length=2) range of axis A-values
#' @param cexPt (numeric)  expansion factor \code{cex} for points (see also \code{\link[graphics]{par}})
#' @param cexSub (numeric)  expansion factor \code{cex} for subtitle (see also \code{\link[graphics]{par}})
#' @param useMar (numeric,length=4) custom margings (see also \code{\link[graphics]{par}})  
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages for debugging 
#' @param silent (logical) suppress messages
#' @return MA-plot only
#' @seealso (for PCA) \code{\link{plotPCAw}}
#' @examples
#' library(wrMisc)
#' set.seed(2005); mat <- matrix(round(runif(600),1),ncol=6)
#' rownames(mat) <- c(rep(letters[1:25],each=3),letters[2:26])
#' MAplotW(mat[,2]-mat[,1], rowMeans(mat))
#' ## assume 2 groups with 3 samples each
#' matMeans <- rowGrpMeans(mat, gr=gl(2,3,labels=LETTERS[3:4]))
#' MAplotW(matMeans[,2]-matMeans[,1], matMeans) 
#' ## assume 2 groups with 3 samples each and run moderated t-test (from package 'limma')
#' tRes <- moderTest2grp(mat,gl(2,3))
#' MAplotW(tRes$Mval, tRes$Amean)                          
#' MAplotW(M=tRes$Mval, A=tRes$means, FCth=1.3) 
#' MAplotW(tRes)
#' MAplotW(tRes, limM=c(-2,2), FCth=1.3) 
#' 
#' @export
MAplotW <- function(Mvalue, Avalue=NULL, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, subTxt=NULL, grayIncrem=TRUE,
  compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=0.7, limM=NULL, limA=NULL, cexPt=NULL, cexSub=NULL, useMar=c(6.2,4,4,2), callFrom="", silent=FALSE,debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="MAplotW")
  opar <- graphics::par(no.readonly = TRUE) 
  on.exit(graphics::par(opar$mar)) 
  on.exit(graphics::par(opar$cex.main)) 
  on.exit(graphics::par(opar$cex.lab)) 
  on.exit(graphics::par(opar$las)) 
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(Avalue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)        # grey
  fcRGB <- c(1,0,0)               # red        for points passing  FC filt line
  if(debug) silent <- FALSE
  if(length(Mvalue) <1) message(" nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if("MArrayLM" %in% class(Mvalue)) {
      ## look for M values
      Mcol <- wrMisc::naOmit(match(c("mvalue","mval","m"), tolower(names(Mvalue))))
      if(length(Mcol) <1) {
        if("means" %in% names(Mvalue)) {
          Mvalue$Mval <- Mvalue$means[,1] - Mvalue$means[,2]
          Mcol <- which(names(Mvalue)=="Mval")
          if(!silent & ncol(Mvalue$means) >2) message("Could not find explicit field for M-values, using 1st and 2nd of group-means for calulating M-values")
        } else stop("Could not find suitable field for M-values in '",namesIn[1],"' !") }    
      ## look for A values
      if(length(Avalue) >0) if(length(as.numeric(Avalue)) != length(as.numeric(Mvalue[[Mcol]]))) {
        if(!silent) message(fxNa,"Content of 'Avalue' seems not to fit to 'Mvalue', ignore and try to use A-values from '",namesIn[1],"' ")
        Avalue <- NULL
      }
      if(length(Avalue) <1) {
        ## no explicit Avalue, try to extract from MArrayLM-object (Mvalue)
        Acol <- wrMisc::naOmit(match(c("amean","avalue","aval","a"),tolower(names(Mvalue))))
        if(length(Acol) >0) Avalue <- Mvalue[[Acol[1]]] else {
          if("means" %in% names(Mvalue)) {
            Avalue <- rowMeans(Mvalue$means[,1:2])
            if(!silent & ncol(Mvalue$means) >2) message("Could not find explicit field for A-values,",
              " using 1st and 2nd of group-means for calulating A-values")
            } else stop("Could not find suitable field for A-values in '",namesIn[1],"' !") }
        if(length(dim(Avalue)) >1) { ANa <- rownames(Avalue)
          Avalue <- as.numeric(Avalue)
          if(length(ANa) >0) names(Avalue) <- ANa }
      } 
      ## recuperate filtering - if present, but only when no custom filtering provided
      if(length(filtFin) <1 | identical(filtFin,FALSE)) {
        Fcol <- wrMisc::naOmit(match(c("filtfin","filter","filt","finfilt"),tolower(names(Mvalue))))
        if(length(Fcol) >0) filtFin <- Mvalue[[Fcol[1]]]
      }
      ## recuperate M values (& dismiss rest of MArrayLM-object)  
      Mvalue <- Mvalue[[Mcol[1]]]
      if(length(dim(Mvalue)) >1) {MNa <- rownames(Mvalue)
        Mvalue <- as.numeric(Mvalue)
        if(length(MNa) >0) names(Mvalue) <- MNa }
      chAM <- length(Mvalue)==length(Avalue)  
      if(!chAM & !silent) message("A- and M- values have different length !!  (M=",length(Mvalue)," vs A=",length(Avalue),")")      
    } else {
      ## thus argument 'Mvalue' is not 'MArrayLM'-object
      ## ... case of explicit Avalue argument
      if(length(Avalue) <1) stop(" argument 'Avalue' is required (if 'Mvalue' not 'MArrayLM'-type object) !")
      if(length(dim(Avalue)) >1) if(ncol(Avalue) >1) {
        if(!silent) message(fxNa," Note, ",namesIn[2]," has ",ncol(Avalue)," columns, using overall mean")
        Avalue <- rowMeans(Avalue,na.rm=TRUE) }
    }
    ## check for (same) order, adjust Mvalue & Avalue according to names
    chNa <- list(MNa=if(length(dim(Mvalue)) >1) rownames(Mvalue) else names(Mvalue),
      ANa=if(length(dim(Avalue)) >1) rownames(Avalue) else names(Avalue))
    nIni <- c(M=length(Mvalue),A=length(Avalue))
    if(length(chNa$MNa) >0 & length(chNa$ANa) >0) {    # ie both have names, so one can match names
      if(!identical(chNa$MNa,chNa$ANa)) {
        matchNa <- wrMisc::naOmit(match(chNa$MNa,chNa$ANa))
        if(length(matchNa) <1) stop("Both 'Mvalue' and 'Avalue' have names, but none of them match !!")
        Avalue <- Avalue[matchNa]
        Mvalue <- wrMisc::naOmit(Mvalue[match(names(Avalue),names(Mvalue))])
      } } else {
        if(length(Mvalue) != length(Avalue)) stop("A- and M- values have different length, but no names to match !!  (M=",length(Mvalue)," vs A=",length(Avalue),")")
      }
    merg <- data.frame(ID=NA,Mvalue=Mvalue,Avalue=Avalue)
    if(length(names(Mvalue)) >0) merg[,1] <- names(Mvalue) else if(length(names(Avalue)) >0) merg[,1] <- names(Avalue)
    ## integrate filtering
    if(length(filtFin) >0) {
      ## if filtFin is matrix use each line with min 1 instance of TRUE,
      if(length(dim(filtFin)) >1) filtFin <- rowSums(as.logical(as.matrix(filtFin))) >0    # use rows with >= 1 TRUE
      if(length(names(filtFin)) >0) {
        matchNa <- wrMisc::naOmit(match(rownames(merg),names(filtFin)))
        if(length(matchNa)==nrow(merg)) merg <- cbind(merg,filtFin=filtFin[matchNa])
      } else if(length(filtFin)==nrow(merg)) merg <- cbind(merg,filtFin=filtFin[matchNa]) 
    } else filtFin <- rep(TRUE,nrow(merg)) 
    if(debug) message(fxNa," ++ DONE extracing columns : ",wrMisc::pasteC(colnames(merg),quo="'"))
    ##
    msg <- " data provided in 'Mvalue' and 'Avalue' "
    if(!silent & nrow(merg) < round(length(Mvalue)/10)) message(" .. note : less than 10% of",msg," were matched") else {
      if(!silent & nrow(merg) < length(Mvalue)/2) message(" .. NOTE : less than 50% of",msg," were matched !!")}
    if(debug) message(msg," were matched to ",nrow(merg)," common entries")
    ## apply filtering
    if(length(filtFin) >0 & !identical(filtFin,FALSE)) {   #  allow filtering provided
      if(sum(filtFin) >0 & sum(filtFin) < nrow(merg)) { merg <- merg[which(merg$filtFin),]
        if(!silent) message(fxNa," filtered (based on 'filtFIn') from ",length(filtFin)," to  ",nrow(merg)," lines")}
    } else filtFin <- TRUE
    nIDco <- sum(c("ID","nredID","uniqID") %in% colnames(merg))                   #  number of heading columns in 'merg'
    Mvalue <- as.numeric(if("Mvalue" %in% colnames(merg)) merg[,"Mvalue"] else merg[,nIDco+1])
    Avalue <- as.numeric(if("Avalue" %in% colnames(merg)) merg[,"Avalue"] else {
      if(length(dim(Mvalue)) >0) merg[,ncol(Mvalue)+nIDco+1] else merg[,nIDco+2]})
    if("Lfdr" %in% colnames(merg)) FdrList <- merg[,"Lfdr"] else {
      if("lfdr" %in% colnames(merg)) FdrList <- merg[,"lfdr"]}
    if(is.null(cexSub)) cexSub <- cexLa +0.05  
    xLab <- "M-value (relative change as log2-ratio)"
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","MA-plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "Volcano-plot"
    alph <- sort(c(0.14,round(0.65/log10(length(Mvalue)),2),0.9))[2]        # alph <- round(12/sqrt(nrow(eBayesLst$Avalue)),2)
      alph2 <- max(round(4/(10+sum(filtFin)^0.25),2),alph)                  # 
    useCex <- max(round(0.7+ 1.4/(1+sum(filtFin,na.rm=TRUE))^0.3,2),0.7)
    ## start plotting
    graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa,cex.lab=cexLa,las=1 )
    baseCol <- if(grayIncrem) grDevices::rgb(0.95,0.95,0.95) else grDevices::rgb(0.7,0.7,0.7)
    graphics::plot(Avalue,Mvalue,pch=16,cex=useCex,main=tit1,xlab="A-value (average log2 intensity)",col=baseCol,ylab=xLab,cex.lab=cexLa,ylim=limM,xlim=limA)
    graphics::abline(h=0,lty=2,col="seagreen")
    if(sum(filtFin,na.rm=TRUE) <1){
      if(!silent) message(fxNa," 0 elements/lines passing 'filtFin' !!")
      FCthrs <- NULL }
    if(length(FCthrs) >0) {
      FCthrs <- log2(FCthrs)
      graphics::abline(h=c(-1,1)*FCthrs+c(0.01,-0.02), col=grDevices::rgb(0.87,0.72,0.72),lty=2)}
    if(grayIncrem & sum(filtFin,na.rm=TRUE) >0) { if(length(FCthrs) <1) {       # no FCthrs, grey incement only on filtFin
      graphics::points(Avalue[filtFin],Mvalue[filtFin],pch=16,cex=useCex,col=grDevices::rgb(0.5,0.5,0.5,alph))
    } else  {
      passFC <- which(filtFin & abs(Mvalue) > FCthrs)
      alph2 <- max(round(1.3/(10+sum(filtFin)^0.25),2),alph+0.05)
      if(!silent) message(" n=",length(Mvalue),"  FCthrs=",signif(2^FCthrs,2),"  filt=",sum(filtFin),
        "  passFC=",length(passFC),"  range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=T),3)),"  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
      if(length(passFC) >0) {
        if(length(passFC) < sum(filtFin)) graphics::points(Avalue[-passFC],Mvalue[-passFC],pch=16,cex=useCex,col=grDevices::rgb(0.5,0.5,0.5,alph))
        graphics::points(Avalue[passFC],Mvalue[passFC],col=grDevices::rgb(0.8,0.01,0.01,alph2),pch=16,cex=useCex)        # red point (only passing filtFin)
      } else graphics::points(Avalue[filtFin],Mvalue[filtFin],pch=16,cex=useCex,col=grDevices::rgb(0.5,0.5,0.5,alph))
      graphics::mtext(paste(if(!is.null(subTxt)) paste("plot",subTxt,"data; "),"n =",length(Mvalue),"; ",
        length(passFC),"(red) point(s) passing filtering (FCthr=",signif(2^FCthrs,2),")"),cex=0.75,line=0.2)
    }} }}
    
