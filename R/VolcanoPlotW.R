#' Volcano-plot (Statistical Test Outcome versus Relative Change)   
#'
#' This type of plot is very common in high-throughput biology, see \href{https://en.wikipedia.org/wiki/Volcano_plot_(statistics)}{Volcano-plot}.
#' Basically, this plot allows comparing the outcome of a statistical test to the differential of the group means (ie log fold-change),
#' 
#' In high-throughput biology data are typically already transformed to log2 and thus, the 'M'-value represents a relative change.
#' Besides, output from statistical testing by \code{\link[wrMisc]{moderTest2grp}} or \code{\link[wrMisc]{moderTestXgrp}} can be directly read to produce Volcano plots for diagnostic reasons.
#' Please note, that plotting a very number of points in transparency (eg >10000) may take several seconds.
#'	 
#' @param Mvalue (numeric or matrix) data to plot; M-values are typically calculated as difference of log2-abundance values and 'pValue' the mean of log2-abundance values;
#'   M-values and p-values may be given as 2 columsn of a matrix, in this case the argument \code{pValue} should remain NULL 
#' @param pValue (numeric, list or data.frame) if \code{NULL} it is assumed that 2nd column of 'Mvalue' contains the p-values to be used
#' @param useComp (integer) choice of one of multiple comparisons present in \code{Mvalue} (if generated using \code{moderTestXgrp()})  
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied 
#' @param ProjNa (character) custom title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC), default at 1.5, set to \code{NA} for omitting
#' @param FdrList (numeric) FDR data or name of list-element
#' @param FdrThrs (numeric) FDR threshold (display as line), default at 0.05, set to \code{NA} for omitting 
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points as increased shades of gray
#' @param col (character) custom color(s) for points of plot (see also \code{\link[graphics]{par}})
#' @param pch (integer) type of symbol(s) to plot (default=16) (see also \code{\link[graphics]{par}}) 
#' @param compNa (character) names of groups compared
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace  
#' @param cexMa (numeric) font-size of title, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexLa (numeric) size of axis-labels, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}}) 
#' @param limM (numeric, length=2) range of axis M-values 
#' @param limp (numeric, length=2) range of axis p-values
#' @param annotColumn (character) column names of annotation to be extracted (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param annColor (character or integer) colors for specific groups of annoatation (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param cexPt (numeric) size of points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexSub (numeric) size of subtitle, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexTxLab (numeric) size of text-labels for points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param namesNBest (integer or character) number of best points (by pValue to display names in figure); if 'passThr' all points passing FDR and FC-filtes will be selected
#' @param NbestCol (character or integer) colors for text-labels of best points
#' @param useMar (numeric,length=4) custom margings (see also \code{\link[graphics]{par}})
#' @param returnData (logical) optional returning data.frame with (ID, Mvalue, pValue, FDRvalue, passFilt) 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages for debugging 
#' @return MA-plot only
#' @seealso (for PCA) \code{\link{plotPCAw}}
#' @examples
#' library(wrMisc)
#' set.seed(2005); mat <- matrix(round(runif(900),2), ncol=9)
#' rownames(mat) <- paste0(rep(letters[1:25],each=4), rep(letters[2:26],4))
#' mat[1:50,4:6] <- mat[1:50,4:6] + rep(c(-1,1)*0.1,25)
#' mat[3:7,4:9] <- mat[3:7,4:9] + 0.7
#' mat[11:15,1:6] <- mat[11:15,1:6] - 0.7
#' ## assume 2 groups with 3 samples each
#' gr3 <- gl(3,3,labels=c("C","A","B"))
#' tRes2 <- moderTest2grp(mat[,1:6], gl(2,3), addResults = c("FDR","means"))
#' # Note: due to the small number of lines only FDR chosen to calculate 
#' VolcanoPlotW(tRes2)
#' ## Add names of points passing custom filters
#' VolcanoPlotW(tRes2, FCth=1.3, FdrThrs=0.2, namesNBest="passThr")
#'
#' ## assume 3 groups with 3 samples each
#' tRes <- moderTestXgrp(mat, gr3, addResults = c("FDR","means"))
#' # Note: due to the small number of lines only FDR chosen to calculate 
#' VolcanoPlotW(tRes)
#' VolcanoPlotW(tRes, FCth=1.3, FdrThrs=0.2)
#' VolcanoPlotW(tRes, FCth=1.3, FdrThrs=0.2, useComp=2)
#'  
#' @export
VolcanoPlotW <- function(Mvalue, pValue=NULL, useComp=1, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, FdrList=NULL, FdrThrs=NULL,
  subTxt=NULL, grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limp=NULL,
  annotColumn=c("SpecType","ProteinName","Accession","Species","Contam","Description"), annColor=NULL, cexPt=NULL, cexSub=NULL, 
  cexTxLab=0.7, namesNBest=NULL, NbestCol=1, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE,debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="VolcanoPlotW")
  opar <- graphics::par(no.readonly=TRUE) 
  on.exit(graphics::par(opar$mar)) 
  on.exit(graphics::par(opar$cex.main)) 
  on.exit(graphics::par(opar$cex.lab)) 
  on.exit(graphics::par(opar$las)) 
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(pValue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)           # grey
  fcRGB <- c(1,0,0)                  # red        for points passing  FC filt line
  splNa <- annot <- ptType <- colPass <- ptBg <- grpMeans <- NULL      # initialize  
  if(debug) silent <- FALSE
  if(length(Mvalue) <1) message(" nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if(length(cexTxLab) <0) cexTxLab <- 0.7
    if("MArrayLM" %in% class(Mvalue)) {
      ## no pValues provided : extract from MArrayLM-object (Mvalue)
      if(length(pValue) <1) {
        ## no explicit pValue, try to extract from MArrayLM-object (Mvalue) => take pValue from MArrayLM-object 
        pcol <- wrMisc::naOmit(match(c("p.value","pvalue","pval","p"), tolower(names(Mvalue))))
        if(length(pcol) >0) pValue <- Mvalue[[pcol[1]]] else stop("can't find suitable element for p-Values")
        if(length(dim(pValue)) >0) if(colnames(pValue)[1]=="(Intercept)" & ncol(pValue) >1) {
          ## extract 2nd col if result from wrMisc::moderTest2grp() 
          pNa <- rownames(pValue)
          pValue <- as.numeric(pValue[,2])
          names(pValue) <- pNa 
        } else {
          ## need to find corresponding of multiple comparisons
          names(useComp) <- colnames(pValue)[useComp]
          splNa <- unlist(strsplit(colnames(pValue)[useComp], "-"))
          pNa <- rownames(pValue)
          if(ncol(pValue) >1) { if(useComp > ncol(pValue)) { useComp <- 1
            names(useComp) <- colnames(pValue)[useComp]
            if(!silent) message(fxNa," argument 'useComp' for pValues invalid or too high; reset to 1") }
          pValue <- as.numeric(pValue[,useComp])
          if(length(pNa) >0) names(pValue) <- pNa }}
      }      
      ## look for M-values (need to create if not available - using useComp checked when extracting pValue)
      Melem <- wrMisc::naOmit(match(c("mvalues","mvalue","m"), tolower(names(Mvalue))))   # which list-element
      if(length(Melem) >0) { Mvalue$Mval <- Mvalue[[Melem]]
        if(length(dim(Mvalue$Mval)) >0) if(ncol(Mvalue$Mval) >1) {
          if(length(useComp) <1) {message(fxNa," set missing 'useComp' to 1"); useComp <- 1 }
          if(useComp[1] > ncol(Mvalue$Mval)) { if(!silent) message(fxNa," re-set invalid 'useComp' to 1"); useComp <- 1 }
          Mvalue$Mval <- Mvalue$Mval[,useComp[1]]}
      } else {
        if("means" %in% names(Mvalue)) {
          if(debug) message(fxNa, "Reconstructing Mvalues based on group-means")
          if(length(splNa) <1) splNa <- NULL
          if(ncol(Mvalue$means) >2 & length(splNa) >1) {
            chMe <- lapply(splNa, function(x) colnames(Mvalue$means) %in% x)
            if(any(sapply(chMe,sum,na.rm=TRUE) <1)) stop("Can't find/match columns of means designated by 'useComp' : ",
              wrMisc::pasteC(useComp[which(sapply(chMe,sum,na.rm=TRUE) <1)], quoteC="'"))
            if(any(sapply(chMe,sum,na.rm=TRUE) >1) & !silent) message(fxNa," Some expected column-names of Mvalue$means appear multiple times ! (using first)")
            grpMeans <- cbind(mean1=Mvalue$means[,which(chMe[[1]])[1]], mean2=Mvalue$means[,which(chMe[[2]])[1]])
          } else grpMeans <- Mvalue$means 
          Mvalue$Mval <- grpMeans[,2] - grpMeans[,1]
          Melem <- which(names(Mvalue)=="Mval")
        } else stop("Could not find suitable field '$means' in '",namesIn[1],"' to construct M-values !")    
      }
      ## pValues : if provided separately, check Corresp p-values & Mvalues
      if(length(pValue) >0) if(length(as.numeric(pValue)) != length(as.numeric(Mvalue[[Melem]]))) {
        if(!silent) message(fxNa,"Content of 'pValue' seems not to fit to 'Mvalue', ignore and try to use p-values from '",namesIn[1],"' ")
        pValue <- NULL
      }
      ## extract FDR from MArrayLM-object 
      if(length(FdrList) <1) {
        ## no explicit pValue, try to extract from MArrayLM-object (Mvalue)
        Fcol <- wrMisc::naOmit(match(c("FDR","BH","lfdr","BY"), tolower(names(Mvalue))))
        if(length(Fcol) >0) { FDRvalue <- Mvalue[[Fcol[1]]]
          ## extract 2nd col if result from wrMisc::moderTest2grp() 
          if(length(dim(FDRvalue)) >0) if(colnames(FDRvalue)[1]=="(Intercept)" & ncol(FDRvalue) >1) {
            pNa <- rownames(FDRvalue)
            FDRvalue <- as.numeric(FDRvalue[,2])
            names(FDRvalue) <- pNa }
        } else {
          FDRvalue <- stats::p.adjust(pValue) 
          if(!silent) message(fxNa,"No FDR data found, generating BH-FDR")}
        ## need to find corresponding of multiple comparisons
        if(length(dim(FDRvalue)) >0) {
          pNa <- rownames(FDRvalue)
          if(ncol(FDRvalue) >1) { if(useComp > ncol(FDRvalue)) { useComp <- 1
            if(!silent) message(fxNa," argument 'useComp' for FDRvalues invalid or too high; reset to 1") }
          FDRvalue <- as.numeric(FDRvalue[,useComp])
          if(length(pNa) >0) names(FDRvalue) <- pNa }}
      }
      ## recuperate filtering - if present, but only when no custom filtering provided
      if(length(filtFin) <1 | identical(filtFin, FALSE)) {
        Fcol <- wrMisc::naOmit(match(c("filtfin","filter","filt","finfilt"),tolower(names(Mvalue))))
        if(length(Fcol) >0) filtFin <- Mvalue[[Fcol[1]]]
      }
      ## recuperate $annot if present and use for symbol
      if("annot" %in% names(Mvalue)) {
        useAnnCol <- match(annotColumn,colnames(Mvalue$annot))      
        if(!is.na(useAnnCol[1])) {                     # annotation (for multiple groups) exists
          ptType <- Mvalue$annot[,useAnnCol[1]]
          if(length(pch) < length(pValue) & length(unique(ptType)) >1) {
            if(length(pch) >1 & !silent) message(fxNa," (invalid pch) using default 'pch' oriented by $annot and starting from 15")
            pch <- 14 + as.numeric(as.factor(ptType))
          } 
          useAnnCol <- wrMisc::naOmit(useAnnCol)
          annot <- Mvalue$annot[,useAnnCol] 
      } }
      if(length(pch)==1) pch <- rep(as.integer(pch), length(pValue))        
      
      ## recuperate M values (& dismiss rest of MArrayLM-object)  
      Mvalue <- Mvalue[[Melem[1]]]
      if(length(dim(Mvalue)) >1) {MNa <- rownames(Mvalue)
        Mvalue <- as.numeric(Mvalue)
        if(length(MNa) >0) names(Mvalue) <- MNa }
      chpM <- length(Mvalue)==length(pValue)  
      if(!chpM & !silent) message("p- and M- values have different length !!  (M=",length(Mvalue)," vs p=",length(pValue),")")
      ## done with extracing MArrayLM-object
      if(!silent) message(fxNa,"Successfully extracted  ",length(Mvalue)," Mvalues and  ",length(pValue)," pValues", if(length(annot) >0) c(" plus anotation"))      
    } else {
      ## thus argument 'Mvalue' is not 'MArrayLM'-object
      ## ... case of explicit pValue argument
      if(length(pValue) <1) stop(" argument 'pValue' is required (if 'Mvalue' not 'MArrayLM'-type object) !")
      if(length(dim(pValue)) >1) if(ncol(pValue) >1) {
        if(!silent) message(fxNa," Note, ",namesIn[2]," has ",ncol(pValue)," columns, using last column")
        pNa <- rownames(pValue)
        pValue <- as.numeric(pValue[,ncol(pValue)] )  
        names(pValue) <- pNa} 
      if(length(FdrList) <1) {
        FDRvalue <- NULL
      }  
    } 
    
    ## need to introduce -log10 to pValue
    chNA <- is.na(pValue)
    if(all(chNA)) stop(fxNa," All p-values are NA, nothing to draw !")
    pValue <- -log10(pValue) 
    ## check for (same) order, adjust Mvalue & pValue according to names
    chNa <- list(MNa=if(length(dim(Mvalue)) >1) rownames(Mvalue) else names(Mvalue),
      pNa=if(length(dim(pValue)) >1) rownames(pValue) else names(pValue))
    nIni <- c(M=length(Mvalue),p=length(pValue))
    if(length(chNa$MNa) >0 & length(chNa$pNa) >0) {        # ie both have names, so one can match names
      if(!identical(chNa$MNa,chNa$pNa)) {
        matchNa <- wrMisc::naOmit(match(chNa$MNa,chNa$pNa))
        if(length(matchNa) <1) stop("Both 'Mvalue' and 'pValue' have names, but none of them match !!")
        pValue <- pValue[matchNa]
        Mvalue <- wrMisc::naOmit(Mvalue[match(names(pValue),names(Mvalue))])
      } } else {
        if(length(Mvalue) != length(pValue)) stop("p- and M- values have different length, but no names to match !!  (M=",length(Mvalue)," vs p=",length(pValue),")")
      }
    if(length(grpMeans) <1) grpMeans <- matrix(rep(NA,2*length(Mvalue)), ncol=2, dimnames=list(names(Mvalue),c("mean1","mean2")))
    merg <- if(length(annot) >0) data.frame(ID=NA, grpMeans, Mvalue=Mvalue, pValue=pValue, FDR=if(length(FDRvalue) >0) FDRvalue else rep(NA,length(pValue)), 
      filtFin=rep(TRUE,length(pValue)), annot, stringsAsFactors=FALSE) else {data.frame(ID=NA, grpMeans, Mvalue=Mvalue, pValue=pValue, stringsAsFactors=FALSE) }
    if(length(names(Mvalue)) >0) merg[,1] <- names(Mvalue) else {if(length(names(pValue)) >0) merg[,1] <- names(pValue)}

    ## adjust col & pch
    if(!any(c(1,length(Mvalue)) %in% length(pch))) {
      if(!silent) message(fxNa,"argument 'pch' should be either length=1 or correspond to length of data, reset to default=16")
      pch <- 16 }
    if(length(col) >1 & length(col) <length(Mvalue)) {
      if(!silent) message(fxNa,"argument 'col' should be either length=1 or correspond to length of data, reset to default=NULL")
      col <- NULL }
    ## integrate filtering
    if(length(filtFin) >0) {
      ## if filtFin is matrix use each line with min 1 instance of TRUE,
      if(length(dim(filtFin)) >1) filtFin <- rowSums(as.logical(as.matrix(filtFin)[,useComp])) >0    # use rows with >= 1 TRUE
      if(length(names(filtFin)) >0) {
        matchNa <- wrMisc::naOmit(match(rownames(merg), names(filtFin)))       
        if(length(matchNa)==nrow(merg)) merg[,"filtFin"] <- filtFin[matchNa]
      } else if(length(filtFin)==nrow(merg)) merg[,"filtFin"] <- filtFin        # no proof that order of filtFin is correct
    } else filtFin <- rep(TRUE, nrow(merg)) 
    if(debug) message(fxNa," ++ DONE extracing columns : ",wrMisc::pasteC(colnames(merg),quo="'"))
    ##
    msg <- " data provided in 'Mvalue' and 'pValue' "
    if(!silent & nrow(merg) < round(length(Mvalue)/10)) message(" .. note : less than 10% of",msg," were matched") else {
      if(!silent & nrow(merg) < length(Mvalue)/2) message(" .. NOTE : less than 50% of",msg," were matched !!")}
    if(debug) message(msg," were matched to ",nrow(merg)," common entries")
    
    ## apply filtering
    if(length(filtFin) >0 & !identical(filtFin,FALSE)) {                         #  use filtering provided
      if(sum(filtFin) >0 & sum(filtFin) < nrow(merg)) { 
        whFilt <- which(merg$filtFin)
        if(length(pch) >1) pch <- pch[whFilt]
        if(length(col) >1) col <- col[whFilt]
        merg <- merg[whFilt,]
        if(!silent) message(fxNa," filtered (based on 'filtFin') from ",length(filtFin)," to  ",nrow(merg)," lines")
        }
    } else filtFin <- rep(TRUE,nrow(merg))
    nIDco <- sum(c("ID","nredID","uniqID") %in% colnames(merg))                   #  number of heading columns in 'merg'
    Mvalue <- as.numeric(if("Mvalue" %in% colnames(merg)) merg[,"Mvalue"] else merg[,nIDco+1])
    pValue <- as.numeric(if("pValue" %in% colnames(merg)) merg[,"pValue"] else {
      if(length(dim(Mvalue)) >0) merg[,ncol(Mvalue) +nIDco +1] else merg[,nIDco+2]})
    if("Lfdr" %in% colnames(merg)) FdrList <- merg[,"Lfdr"] else {
      if("lfdr" %in% colnames(merg)) FdrList <- merg[,"lfdr"]}

    ## prepare for  plotting
    if(is.null(cexSub)) cexSub <- cexLa +0.05  
    xLab <- "M-value (log2 fold-change)"
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","Volcano-plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "Volcano-plot"
    if(length(FCthrs) <1) FCthrs <- 1.5 
    if(length(FdrThrs) <1) FdrThrs <- 0.05 
    ## count no of passing
    passFC <- if(length(FCthrs) ==1 & !any(is.na(FCthrs))) abs(Mvalue) > log2(FCthrs) else filtFin      ## convert FCthrs to log2
    passFdr <- if(length(FdrThrs) ==1 & !any(is.na(FdrThrs))) FDRvalue <= FdrThrs else filtFin
    passAll <- filtFin & passFC & passFdr
    chNA <- is.na(passAll)                              # passFdr may contain NAs
    if(any(chNA)) passAll[which(chNA)] <- FALSE 
    if(debug) message(fxNa,"  ",sum(passFC,na.rm=TRUE)," passing FCthrs ; ",sum(passFdr,na.rm=TRUE)," passing FdrThrs ; combined ",sum(passAll,na.rm=TRUE))
    ## color for points passing filter
    if(length(col) <1) {
      alph <- sort(c(0.14, round(0.6/log10(length(Mvalue)),2), 0.8))[2]       # alph <- round(12/sqrt(nrow(eBayesLst$pValue)),2)
          alph2 <- max(round(4/(10 +sum(filtFin)^0.25),2), alph) 
          alph3 <- max(round(1.3/(10 +sum(filtFin)^0.25),2), alph +0.05)                 # alternative for alph2 
      useCol <- if(grayIncrem) grDevices::rgb(0.35,0.35,0.35,alph) else grDevices::rgb(0.7,0.7,0.7)
      useCex <- if(length(cexPt) >0) cexPt else max(round(0.8 +2/(1 +sum(filtFin,na.rm=TRUE))^0.28,2), 1.1)
      colPass <- grDevices::rgb(0.8,0.01,0.01, alph2)                 # (default) red
      whCol <- rep(1,length(Mvalue))
      if(any(passAll)) whCol[which(passAll)] <- 2                     # assign color for those passing
      useCol <- c(useCol,colPass)[whCol]  
    } else useCol <- col
    ## adjust fill color for open symbols
    chPch <- pch %in% c(21:25)
    if(any(chPch)) {ptBg <- col ; col[which(chPch)] <- grDevices::rgb(0.2,0.2,0.2, max(alph,0.4)) }
    
    ## main graphic
    #graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa, cex.lab=cexLa, las=1)
    graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa, las=1)
    #graphics::plot(Mvalue, pValue, pch=pch, cex=useCex, main=tit1, ylab="- log10 p-value (uncorrected)", col=useCol,xlab=xLab,cex.lab=cexLa,xlim=limM,ylim=limp)
    graphics::plot(Mvalue, pValue, pch=pch, cex=useCex, main=tit1, ylab="- log10 p-value (uncorrected)", col=useCol, xlab=xLab,cex.lab=cexLa,xlim=limM,ylim=limp,pt.bg=ptBg)
    sTxt <- paste(if(length(splNa) >0) paste(names(useComp),": "), if(!is.null(subTxt)) paste("plot",subTxt,"data; "), "n =",length(Mvalue),
      if(!all(is.na(c(FCthrs,FdrThrs)))) paste("; ",sum(passAll, na.rm=TRUE),"(red) points passing thresholds (FCthr=",as.character(FCthrs),
        ", FdrThrs=",as.character(FdrThrs),")"))
    graphics::mtext(sTxt,cex=0.75,line=0.2)
    if(!all(is.na(c(FCthrs,FdrThrs)))) { 
      if(debug) message(fxNa," n=",length(Mvalue),"  FCthrs=",as.character(FCthrs),"  filt.ini=", sum(filtFin, na.rm=TRUE),
        "  passAll=",sum(passAll,na.rm=TRUE)," ; range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=TRUE),3))," ;  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
      graphics::abline(v=c(-1,1)*log2(FCthrs) +c(0.01,-0.02), col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
    if(!any(is.na(FdrThrs))) {
      dThr <- FDRvalue -FdrThrs
      dThr2 <- which(dThr <0) 
      if(length(dThr2) >0) {
        limP <- dThr2[which(dThr[dThr2] ==max(dThr[dThr2], na.rm=TRUE))]
        limP <- min(pValue[limP], na.rm=TRUE)
        graphics::abline(h=limP -0.015, col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
    }

    ## add names to best points
    if(length(namesNBest) >0) { 
      if(identical(namesNBest,"passThr") | identical(namesNBest,"signif")) namesNBest <- sum(passAll) 
      if(!is.integer(namesNBest)) namesNBest <- try(as.integer(namesNBest))
      tmP <- as.numeric(merg[,"pValue"])
      if(any(!passAll)) tmP[which(!passAll)] <- min(tmP,na.rm=TRUE) -1
      useLi <- order(tmP, decreasing=TRUE)[1:namesNBest]      
      xOffs <- signif(abs(diff(if(length(limM)==2) limp else range(Mvalue, na.rm=TRUE) ))/130,3)
      yOffs <- signif(abs(diff(if(length(limp)==2) limp else range(pValue, na.rm=TRUE) ))/90,3)
      ## alternative names to display
      dispNa <- merg[useLi,1] 
      if(length(annotColumn) >2) if(annotColumn[2] %in% colnames(annot)) {          # if available, use 2nd annotColumn for text lables, default annot-column 'ProteinName'
        dispNa <- merg[useLi,annotColumn[2]] }
      if(any(is.na(dispNa))) message(fxNa," Some names for points to highlight are and won't get displayed !")  
      if(length(NbestCol) <1) NbestCol <- 1
      graphics::text(Mvalue[useLi] +xOffs, pValue[useLi] +yOffs, dispNa, cex=cexTxLab, col=NbestCol, adj=0)
    }          

    ## custom colors for significant points
    if(length(annColor) >0 & sum(passAll) >0) {  ## replot points passing thresholds according to colors given
      useLi <- which(passAll)
      if(length(annColor) ==length(pValue)) useCol <- annColor[useLi] else {
        useCol <- annColor[as.numeric(as.factor(merg[,annotColumn[1]]))][useLi]
      } 
      graphics::points(Mvalue[useLi], pValue[useLi], pch=pch[useLi], cex=useCex, col=useCol, pt.bg=ptBg)    
    } else annColor <- NULL
    
    ## legend (if multiple symbols)
    ch1 <- unique(pch)
    if(length(ch1) >1) {
      legLab <- wrMisc::naOmit(if(length(ptType) >0) unique(ptType) )
      legLoc <- wrGraph::checkForLegLoc(cbind(Mvalue, pValue), sampleGrp=legLab, showLegend=FALSE)
      ##
      legCex <- stats::median(c(useCex,cexTxLab,1.2),na.rm=TRUE)
      ptCol <- if(length(colPass) >0) {if(length(annColor) >0 & length(annotColumn) >0) annColor[order(unique(wrMisc::naOmit(merg[,annotColumn[1]])))] else colPass} else useCol[1]
      graphics::legend(legLoc$loc, legend=rev(legLab), col=rev(ptCol), text.col=1, pch=rev(wrMisc::naOmit(unique(pch))), if(length(ptBg) >0) pt.bg=ptBg, 
        cex=legCex, pt.cex=1.2*legCex, xjust=0.5, yjust=0.5)        # as points
    }

  ## export results
  if(returnData) cbind(merg[,1:5],FDRvalue=FDRvalue,merg[,-(1:5)]) 
  } }
    
