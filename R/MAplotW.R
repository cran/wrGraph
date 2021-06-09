#' MA-plot (differential intensity versus average intensity)   
#'
#' This type of plot for display of relative changes versus (mean) absolute abundance is very common in high-throughput biology, see \href{https://en.wikipedia.org/wiki/MA_plot}{MA-plot}.
#' Basically one compares two independent series of measures (ie gene transcript or protein abundance values) of 2 samples/data-sets or the means of 2 groups of replicates.
#' And the log-fold-change ('Minus'=M) is plotted againts the absolute mean value ('Average'=A).
#' Furthermore, output from statistical testing by \code{\link[wrMisc]{moderTest2grp}} or \code{\link[wrMisc]{moderTestXgrp}} can be directly read to produce MA plots for diagnostic purpose.
#' Please note, that plotting a very high number of points in transparency (eg >10000) may take several seconds.
#'	 
#' @param Mvalue (numeric, list or MArrayLM-object) main data to plot; if numeric, the content will be used as M-values (and A-values must be provided separateley);
#'  if list or MArrayLM-object, it must conatin list-elements named \code{Mvalue} and \code{means} to extract all information needed for plotting
#' @param Avalue (numeric, list or data.frame) if \code{NULL} it is assumed that M-values can be extracted form argument \code{Avalue}
#' @param useComp (integer) choice of one of multiple comparisons present in \code{Mvalue} (if generated using \code{moderTestXgrp()})  
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied 
#' @param ProjNa (character) custom title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC)
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points (not exceeding threshold) as increased shades of gray 
#' @param col (character) custom color(s) for points of plot (see also \code{\link[graphics]{par}})
#' @param pch (integer) type of symbol(s) to plot (default=16) (see also \code{\link[graphics]{par}}) 
#' @param compNa depreciated, please use \code{useComp} instead
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace  
#' @param cexMa (numeric) font-size of title, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexLa (numeric) size of axis-labels, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}}) 
#' @param limM (numeric, length=2) range of axis M-values 
#' @param limA (numeric, length=2) range of axis A-values
#' @param annotColumn (character) column names of annotation to be extracted (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot).
#'   The first entry (typically 'SpecType') is used for different symbols in figure, the second (typically 'GeneName') is used as prefered text for annotating the best points (if \code{namesNBest} allows to do so.)
#' @param annColor (character or integer) colors for specific groups of annotation (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param cexPt (numeric) size of points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexSub (numeric) size of subtitle, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexTxLab (numeric) size of text-labels for points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param namesNBest (integer or character, length=1) number of best points to add names in figure; if 'passThr' all points passing FC-filter will be selected; 
#'   if the initial object \code{Mvalue} contains a list-element called 'annot' the second of the column specified in argument \code{annotColumn} will be used as text
#' @param NbestCol (character or integer) colors for text-labels of best points
#' @param NaSpecTypeAsContam (logical) consider lines/proteins with \code{NA} in Mvalue$annot[,"SpecType"] as contaminants (if a 'SpecType' for contaminants already exits)
#' @param useMar (numeric,length=4) custom margings (see also \code{\link[graphics]{par}})
#' @param returnData (logical) optional returning data.frame with (ID, Mvalue, Avalue, FDRvalue, passFilt) 
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging 
#' @return MA-plot only
#' @seealso (for PCA) \code{\link{plotPCAw}}
#' @examples
#' library(wrMisc)
#' set.seed(2005); mat <- matrix(round(runif(600),2), ncol=6)
#' rownames(mat) <- c(rep(letters[1:25],each=3), letters[2:26])
#' MAplotW(mat[,2] -mat[,1], A=rowMeans(mat))
#' ## assume 2 groups with 3 samples each
#' matMeans <- rowGrpMeans(mat, gr=gl(2,3,labels=LETTERS[3:4]))
#' MAplotW(M=matMeans[,2] -matMeans[,1], A=matMeans) 
#' ## assume 2 groups with 3 samples each and run moderated t-test (from package 'limma')
#' tRes <- moderTest2grp(mat, gl(2,3))
#' MAplotW(tRes$Mval, tRes$Amean)                          
#' MAplotW(M=tRes$Mval, A=tRes$means, FCth=1.3) 
#' MAplotW(tRes)
#' MAplotW(tRes, limM=c(-2,2), FCth=1.3) 
#' 
#' @export
MAplotW <- function(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, subTxt=NULL,
  grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limA=NULL,
  annotColumn=c("SpecType","GeneName","EntryName","Accession","Species","Contam"), annColor=NULL, cexPt=NULL, cexSub=NULL, cexTxLab=0.7, 
  namesNBest=NULL, NbestCol=1, NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE,debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="MAplotW")
  opar <- graphics::par(no.readonly=TRUE) 
  on.exit(graphics::par(opar$mar)) 
  on.exit(graphics::par(opar$cex.main)) 
  on.exit(graphics::par(opar$cex.lab)) 
  on.exit(graphics::par(opar$las)) 
  plotByFdr <- TRUE  #
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(Avalue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)           # grey
  fcRGB <- c(1,0,0)                  # red        for points passing  FC filt line
  multiComp <- TRUE                    # initialize  
  splNa <- annot <- ptType <- colPass <- ptBg <- grpMeans <- pcol <- FDRvalue <- FdrList <- FDRty <- NULL      # initialize
  if(length(pch) <1) pch <- 16
  if(debug) silent <- FALSE
  if(debug) message(" length Mvalue ",length(Mvalue)," ; length Avalue ",length(Avalue)," ; useComp ",useComp)
  if(identical(col,"FDR")) {FDR4color <- TRUE; col <- NULL} else FDR4color <- FALSE
  if(length(Mvalue) <1) message(" nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if(length(cexTxLab) <0) cexTxLab <- 0.7
    if("MArrayLM" %in% class(Mvalue) | "list" %in% class(Mvalue)) {
      ## try working based on MArrayLM-object (Mvalue)
      if(debug) message("'",namesIn[1],"' is list or MArrayLM-object ")
      ## initial check of useComp
      if(length(useComp) >1) { useComp <- wrMisc::naOmit(useComp)[1]
        if(!silent) message(fxNa," argument 'useComp' should be integer of length=1; using only 1st entry") }
      if(length(useComp) <1) { useComp <- 1
        if(!silent) message(fxNa," argument 'useComp' invalid, setting to 1") }
      ## address multiple questions (via useComp): need to check 'useComp', thus need to locate FDRvalues or pValues for group-names ..
      pcol <- wrMisc::naOmit(match(c("p.value","pvalue","pval","p"), tolower(names(Mvalue))))          
      FDRcol <- wrMisc::naOmit(match(c("fdr","bh","lfdr","by","bonferroni"), tolower(names(Mvalue))))
      ## extract FDR-values (or p-values if no FDR available)
      if(length(FDRcol) >0) FDRvalue <- Mvalue[[FDRcol[1]]] else if(length(pcol) >0) FDRvalue <- Mvalue[[pcol[1]]] 
      if(FDR4color | useComp >1) {
        if(debug) message(" FDR4color=",FDR4color,"  need to extract FDR or p-values")
        if(length(dim(FDRvalue)) >0) {    # FDRvalues are present as multi-column
          if(colnames(FDRvalue)[1]=="(Intercept)" & ncol(FDRvalue) >1) {
            ## extract 2nd col if result from wrMisc::moderTest2grp()
            if(debug) message(" extract 2nd col if result from wrMisc::moderTest2grp()")
            pNa <- rownames(FDRvalue)
            FDRvalue <- as.numeric(FDRvalue[,2])
            names(FDRvalue) <- pNa
            multiComp <- FALSE 
          } else {
            ## select corresponding of multiple comparisons
            if(debug) message(" select corresponding FDR of multiple comparisons")
            if(useComp > ncol(FDRvalue)) { useComp <- 1
              if(!silent) message(fxNa," argument 'useComp' for FDR-Values invalid or too high; reset to 1") }
            names(useComp) <- colnames(FDRvalue)[useComp]
            pNa <- rownames(FDRvalue)
            FDRvalue <- as.numeric(FDRvalue[,useComp])
            if(length(pNa) >0) names(FDRvalue) <- pNa }
          if(length(FDRcol) <1) { FDRty <- "BH"
            FDRvalue <- stats::p.adjust(FDRvalue, method="BH")   # transform p.value to BH-FDR
          } else FDRty <- names(Mvalue)[FDRcol[1]]
        } else {useComp <- 1}
      } else {     # no need to use useComp 
        plotByFdr <- FALSE
        if(useComp >1) { useComp <- 1 
          if(!silent) message(fxNa," ignoring value of argument 'useComp' since no p-Values or FDR-values found")}
      }
      ## look for M-values (need to create if not available - using useComp checked)
      Melem <- wrMisc::naOmit(match(c("mvalues","mvalue","m"), tolower(names(Mvalue))))      # which list-element

      ## look for group-means & identify column association to current question/pairwise comparison
      if("means" %in% names(Mvalue)) {
        ## identify sample-groups to comparsison(s) - needed lateron
        if(debug) message(" identify sample-groups to comparsison(s)")
        pairwCol <- wrMisc::sampNoDeMArrayLM(Mvalue, useComp, lstMeans="means", lstP=names(Mvalue)[FDRcol[1]], callFrom=fxNa,silent=silent) 
        grpMeans <- cbind(mean1=Mvalue$means[,pairwCol[1]], mean2=Mvalue$means[,pairwCol[2]])  
        ## are all group-means needed (for exporting) ??
      } else warning("Could not find suitable field '$means' in '",namesIn[1],"'")    
      
      if(length(Melem) >0) {            ## M-values are available
        if(debug) message(" M-values are available")
        Mvalue$Mval <- Mvalue[[Melem]]
        if(length(dim(Mvalue$Mval)) >0) if(ncol(Mvalue$Mval) >1) {          
          if(useComp[1] > ncol(Mvalue$Mval)) { if(!silent) message(fxNa," 'useComp' is too high, set to 1"); useComp <- 1 }
          Mvalue$Mval <- Mvalue$Mval[,useComp]}
      } else {                               # need to construct M-values based on means
        if("means" %in% names(Mvalue)) {
          ## construct Mvalue based on means (only one/current pairwise comparison needed)
          if(debug) message(" construct Mvalue based on means")
          Mvalue$Mval <- grpMeans[,2] - grpMeans[,1]                           
          Melem <- which(names(Mvalue)=="Mval")                # update
        } else stop("Can't construct M-values since suitable field '$means' missing in '",namesIn[1],"' !")    
      }

      ## now one can check if 'Avalue' & Mvalue match (only if pValues or FDR values used for coloring)
      if(length(Avalue) >0) if(length(Avalue) != length(Mvalue$Mval)) { Avalue <- NULL
        if(!silent) message(fxNa,"Invalid entry of 'Avalue' (length=",length(Avalue),"but expecting",length(Mvalue$Mval),")")}
      if(length(Avalue) != length(Mvalue$Mval)) Avalue <- Mvalue$Aval <- rowMeans(grpMeans[,1:2], na.rm=TRUE)    # define fresh A-values 

      ## recuperate filtering - if present, but only when no custom filtering provided
      if(length(filtFin) <1 | identical(filtFin, FALSE)) {
        Filcol <- wrMisc::naOmit(match(c("filtfin","filter","filt","finfilt"), tolower(names(Mvalue))))
        filtFin <- if(length(Filcol) >0) Mvalue[[Filcol[1]]] else rep(TRUE,length(Avalue))
        if(length(dim(filtFin)) >1) filtFin <- filtFin[,useComp]
        if(debug) message(" recuperate filtering - if present (length ",length(filtFin),")")
      }

      ## recuperate $annot if present and use for symbol
      if("annot" %in% names(Mvalue)) {
        useAnnCol <- match(annotColumn, colnames(Mvalue$annot))      
        if(!is.na(useAnnCol[1])) {                         # annotation (for multiple groups) exists
          ptType <- Mvalue$annot[,useAnnCol[1]]            # SpecType
          chNA <- is.na(ptType)
          ## associate NAs from 'SpecType' in ptType with conta ?
          if(NaSpecTypeAsContam) {
            chConta <- tolower(ptType) %in% c("contaminant","contam","conta","cont")
            if(any(chConta)) ptType[which(is.na(ptType))] <- unique(ptType[which(chConta)])[1]}          
          if(any(is.na(ptType))) ptType[which(chNA)] <- "NA" 
          if(length(pch) < length(Avalue) & length(unique(wrMisc::naOmit(ptType))) >1) {
            if(length(pch) >1 & !silent) message(fxNa," (invalid pch) using default 'pch' oriented by $annot and starting from 15")
            pch <- 14 + as.integer(as.factor(ptType))            
          } 
          useAnnCol <- wrMisc::naOmit(useAnnCol) }
        annot <- Mvalue$annot[,useAnnCol] 
        if(annotColumn[1] %in% colnames(annot)) annot[,annotColumn[1]] <- ptType   # update with NAs trasformed to "NA"
        }          

      if(length(pch)==1) pch <- rep(as.integer(pch), length(Avalue))
      if(length(filtFin) <1) filtFin <- rep(TRUE, length(Avalue))        
      
      
      ## recuperate M values (& dismiss rest of MArrayLM-object)  
      Mvalue <- Mvalue$Mval
      if(length(dim(Mvalue)) >1) { MNa <- rownames(Mvalue)
        Mvalue <- as.numeric(Mvalue)
        if(length(MNa) >0) names(Mvalue) <- MNa }
      ## additional check for length 
      chpM <- length(Mvalue)==length(Avalue)  
      if(!chpM & !silent) message(fxNa,"trouble ahead ? p- and M- values have different length !!  (M=",length(Mvalue)," vs A=",length(Avalue),")")

      ## done with extracting MArrayLM-object    
      if(!silent) message(fxNa,"Successfully extracted  ",length(Mvalue)," Mvalues and  ",length(Avalue)," Avalues", if(length(annot) >0) c(" plus annotation"))      
    } else {
      ## thus argument 'Mvalue' is not 'MArrayLM'-object
      ## ... case of explicit Avalue argument
      if(length(Avalue) <1) stop(" argument 'Avalue' is required (if 'Mvalue' not 'MArrayLM'-type object) !")
      multiComp <- FALSE
      if(length(dim(Avalue)) >1) if(ncol(Avalue) >1) {
        if(!silent) message(fxNa," Note, ",namesIn[2]," has ",ncol(Avalue)," columns, using last column")
        pNa <- rownames(Avalue)
        Avalue <- as.numeric(Avalue[,ncol(Avalue)] )  
        names(Avalue) <- pNa} 
    }
    
    chNA <- is.na(Avalue)
    if(all(chNA)) stop(fxNa," All A-values are NA, nothing to draw !")

    ## check for (same) order, adjust Mvalue & Avalue according to names
    chNa <- list(MNa=if(length(dim(Mvalue)) >1) rownames(Mvalue) else names(Mvalue),
      ANa=if(length(dim(Avalue)) >1) rownames(Avalue) else names(Avalue))
    nIni <- c(M=length(Mvalue), A=length(Avalue))
    if(length(chNa$MNa) >0 & length(chNa$ANa) >0) {        # ie both have names, so one can match names
      if(!identical(chNa$MNa,chNa$ANa)) {
        matchNa <- wrMisc::naOmit(match(chNa$MNa,chNa$ANa))
        if(length(matchNa) <1) stop("Both 'Mvalue' and 'Avalue' have names, but none of them match !!")
        Avalue <- Avalue[matchNa]
        Mvalue <- wrMisc::naOmit(Mvalue[match(names(Avalue),names(Mvalue))])
      } } else {
        if(length(Mvalue) != length(Avalue)) stop("A- and M- values have different length, but no names to match !!  (M=",length(Mvalue)," vs A=",length(Avalue),")")
      }
    if(length(grpMeans) <1) grpMeans <- matrix(rep(NA,2*length(Mvalue)), ncol=2, dimnames=list(names(Mvalue),c("mean1","mean2")))
    if(length(pch)==1 & length(Mvalue) >1) pch <- rep(pch, length(Mvalue))
    if(length(pch) != length(Mvalue) & (length(pch) >1)) { if(!silent) message(fxNa," bizzare entry for 'pch'")
      pch <- rep(pch, length(Mvalue))[1:length(Mvalue)]}
    ## start creating merged data for plot (& export)
    merg <- if(length(annot) >0) data.frame(ID=NA, grpMeans, Mvalue=Mvalue, Avalue=Avalue, FDR=if(length(FDRvalue) >0) FDRvalue else rep(NA,length(Mvalue)), 
      filtFin=if(length(filtFin) >0) filtFin else rep(TRUE,length(Mvalue)), annot, pch=pch, stringsAsFactors=FALSE) else {
      data.frame(ID=NA, grpMeans, Mvalue=Mvalue, Avalue=Avalue, FDR=if(length(FDRvalue) >0) FDRvalue else rep(NA,length(Mvalue)), filtFin=if(length(filtFin) >0) filtFin else rep(TRUE,length(Mvalue)), pch=pch, stringsAsFactors=FALSE) }
    if(length(names(Mvalue)) >0) merg[,1] <- names(Mvalue) else {if(length(names(Avalue)) >0) merg[,1] <- names(Avalue)}
    ## replace NA in 'SpecType' by 'NA'

    if(annotColumn[1] %in% colnames(merg)) { chNa <- is.na(merg[,annotColumn[1]])     # replace NAs in col "SpecType" by "NA"
      if(any(chNa)) merg[which(chNa),annotColumn[1]] <- "NA"
    } else { merg <- cbind(merg, rep(1,nrow(merg)))            # add colum for 'SpecType'
      colnames(merg)[ncol(merg)] <- annotColumn[1] }

    ## adjust col (color) & pch
    if(!any(c(1,length(Mvalue)) %in% length(pch))) {
      if(!silent) message(fxNa,"argument 'pch' should be either length=1 or correspond to length of data, reset to default=16")
      pch <- 16 }
    if(length(col) >1 & length(col) <length(Mvalue)) {
      if(!silent) message(fxNa,"argument 'col' should be either length=1 or correspond to length of data, reset to default=NULL")
      col <- NULL }

    ## prepare/integrate FILTERING
    if(length(filtFin) >0) {
      ## if filtFin is matrix use each line with min 1 instance of TRUE,
      if(length(dim(filtFin)) >1) filtFin <- as.logical(as.matrix(filtFin)[,useComp])    # use rows with >= 1 TRUE
      if(length(names(filtFin)) >0) {
        matchNa <- wrMisc::naOmit(match(rownames(merg), names(filtFin)))       
        if(length(matchNa)==nrow(merg)) merg[,"filtFin"] <- filtFin[matchNa]
      } else if(length(filtFin)==nrow(merg)) merg[,"filtFin"] <- filtFin        # no proof that order of filtFin is correct
    } else filtFin <- rep(TRUE, nrow(merg)) 
    if(debug) message(fxNa," ++ DONE extracting columns : ",wrMisc::pasteC(colnames(merg),quo="'"))

    ## apply filtering
    msg <- " data provided in 'Mvalue' and 'Avalue' "
    if(!silent & nrow(merg) < round(length(Mvalue)/10)) message(" .. note : less than 10% of",msg," were matched") else {
      if(!silent & nrow(merg) < length(Mvalue)/2) message(" .. NOTE : less than 50% of",msg," were matched !!")}
    if(debug) message(msg," were matched to ",nrow(merg)," common entries")    

    ## apply filtering (keep all lines where at least one condition passes)
    if(length(filtFin) >0 & !identical(filtFin, FALSE)) {                         #  use filtering provided
      if(sum(filtFin) >0 & sum(filtFin) < nrow(merg)) { 
        whFilt <- which(merg$filtFin)
        if(length(pch) >1) pch <- pch[whFilt]
        if(length(col) >1) col <- col[whFilt]
        merg <- merg[whFilt,]
        if(!silent) message(fxNa," filtered (based on 'filtFin') from ",length(filtFin)," to  ",nrow(merg)," lines")
      }
    } else filtFin <- rep(TRUE, nrow(merg))
    
    ## sort merg, so that legend always gets constructed the same order, ascending ('ascend') or descending ('descend')    
    ## update ..
    nIDco <- sum(c("ID","nredID","uniqID") %in% colnames(merg))                   #  number of heading columns in 'merg'
    Mvalue <- as.numeric(if("Mvalue" %in% colnames(merg)) merg[,"Mvalue"] else merg[,nIDco+1])
    Avalue <- as.numeric(if("Avalue" %in% colnames(merg)) merg[,"Avalue"] else {
      if(length(dim(Mvalue)) >0) merg[,ncol(Mvalue) +nIDco +1] else merg[,nIDco+2]})
    if("Lfdr" %in% colnames(merg)) FdrList <- merg[,"Lfdr"] else {
      if("lfdr" %in% colnames(merg)) FdrList <- merg[,"lfdr"]}
    pch <- merg[,"pch"]            # update
    ptType <- if(annotColumn[1] %in% colnames(merg)) merg[,annotColumn[1]] else rep(1,nrow(merg))     # update "SpecType"
    
    ## prepare for  plotting
    if(is.null(cexSub)) cexSub <- cexLa +0.05  
    xLab <- paste("A-value", if(!batchFig) "(average abundance)")
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","MA-plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "MA-plot"
    if(length(FCthrs) <1) FCthrs <- 1.5 
    #needed?#if(length(FdrThrs) <1) FdrThrs <- 0.05 

    ## count no of passing
    passAll <- passFC <- if(length(FCthrs) ==1 & !any(is.na(FCthrs))) abs(merg[,"Mvalue"]) >= log2(FCthrs) else merg[,"filtFin"]      ## convert FCthrs to log2
    passAll <- merg[,"filtFin"] & passFC    

    chNA <- is.na(passAll)                              # passFdr may contain NAs
    if(any(chNA)) passAll[which(chNA)] <- FALSE 
    if(debug) message(fxNa,"  ",sum(passFC,na.rm=TRUE)," passing FCthrs ")
    ## color for points passing filter
    if(length(col) >0) if(length(col) != nrow(merg)) { col <- NULL
      if(!silent) message(fxNa," invalid entry for 'col', should be of length=",nrow(merg),", resetting to default")}

    if(length(col) <1) {
      alph <- sort(c(0.14, round(0.6/log10(length(Mvalue)),2), 0.8))[2]       # alph <- round(12/sqrt(nrow(eBayesLst$pValue)),2)
      alph2 <- sort(c(round(7/(5 +sum(passAll)^0.7),2), alph,0.9))[2]                   # for points passing thresholds
      useCol <- if(grayIncrem) grDevices::rgb(0.35,0.35,0.35,alph) else grDevices::rgb(0.7,0.7,0.7)  # basic color
      useCex <- if(length(cexPt) >0) cexPt else max(round(0.8 + 2/(1 +sum(filtFin, na.rm=TRUE))^0.28,2), 1.1)
      
      chCol <- unique(merg[, annotColumn[1]])      # check how many different colors may be needed
      chNaC <- is.na(chCol)
      if(any(chNaC)) chCol[which(chNaC)] <- "NA"    
      if(length(annColor) >0) {colPass <- annColor} else if(length(chCol) >4) {
        colPass <- cbind(red=c(141,72,90,171, 220,253,244,255), green=c(129,153,194,221, 216,174,109,0), blue=c(194,203,185,164, 83,97,67,0))       
        colPass <- grDevices::rgb(red=colPass[,1], green=colPass[,2], blue=colPass[,3], alph2, maxColorValue=255)
        if(length(chCol) >8) { colPass <- c(colPass, rep(colPass[8], length(chCol) -8))
          if(!silent) message(fxNa," > 8 different groups found, using 8th color after 7th group")}
      } else colPass <- grDevices::rgb(c(0.95,0.2,0,0.75), c(0.15,0.2,0.9,0.35), c(0.15,0.95,0,0.8), alph2)    # red, blue, green, purple (luminosity adjusted) 
      useCol <- rep(useCol[1], nrow(merg))         # fuse basic gray to colors for different types

      ## integrate names of annColor as order of colPass
      if(length(names(annColor)) >0) {
        uniTy <- unique(merg[which(passAll),annotColumn[1]])
        colPass <- colPass[match(names(annColor), uniTy)]
      }     
      ## assign color for those passing
      if(any(passAll)) { if(FDR4color) {
        useCol[which(passAll)] <- .colorByPvalue(merg$FDR[which(passAll)])
        
        } else {
        useCol[which(passAll)] <- colPass[ if(length(unique(merg[which(passAll),annotColumn[1]])) >1) wrMisc::levIndex(
          merg[which(passAll),annotColumn[1]]) else rep(1,sum(passAll))] } }  # assign colors for those passing 
    } else useCol <- col
        
    ## adjust fill color for open symbols
    chPch <- pch %in% c(21:25)
    if(any(chPch)) { ptBg <- useCol
      ptBg[which(chPch)] <- useCol[which(chPch)]    # background color for filled symbols
      useCol[which(chPch)] <- 1                     # contour as black
    }

    ## main graphic
    graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa, las=1)
    ## rather directly plot FDR
    graphics::plot(merg[,"Avalue"], merg[,"Mvalue"], pch=pch, cex=useCex, main=tit1, 
      ylab="M value (log FC)", col=useCol, xlab=xLab, cex.lab=cexLa, xlim=limM,ylim=limA, pt.bg=ptBg)    
    
    sTxt <- if(length(subTxt) ==1) subTxt else { if(multiComp) paste0(if(length(names(useComp)) >0) names(useComp) else paste0("useComp=",useComp),"; ",collapse="")}
    sTxt <- paste0(sTxt,"n=",length(Mvalue),
      if(!all(is.na(c(FCthrs)))) paste(";",sum(passAll, na.rm=TRUE),"(color) points passing", if(!is.na(FCthrs)) paste0("FCthr=", as.character(signif(FCthrs,3))) ))
    graphics::mtext(sTxt, cex=0.75, line=0.2)    
    
    if(!all(is.na(FCthrs))) { 
      if(debug) message(fxNa," n=",length(Mvalue),"  FCthrs=",as.character(FCthrs),"  filt.ini=", sum(filtFin, na.rm=TRUE),"  passAll=",sum(passAll,na.rm=TRUE),
        " ; range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=TRUE),3))," ;  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
      graphics::abline(h=c(-1,1)*(log2(FCthrs) + diff(graphics::par("usr")[1:2])/500), col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
    
    ## add names to best points
    if(length(namesNBest) >0) { 
      if(any(sapply( c("passThr","pass","passFC"), identical, namesNBest))) namesNBest <- sum(passAll)
      if(!is.integer(namesNBest)) namesNBest <- try(as.integer(namesNBest))
      if(namesNBest >0 & any(passAll)) {      
        useLi <- if(any(!passAll)) which(passAll) else 1:nrow(merg)
        tmP <- as.numeric(merg[useLi,"Avalue"])
        names(tmP) <- rownames(merg)[useLi]
        ## look for more informative names to display
        if(length(annot) >0) {
          proNa <- annot[match(names(tmP), rownames(annot)), annotColumn[2]]   # normally 'Description'
          chNa <- is.na(proNa)
          if(!all(chNa)) names(tmP)[which(!chNa)] <- proNa[which(!chNa)]
        }        
        useL2 <- order(tmP, decreasing=TRUE)[1:min(namesNBest,sum(passAll))]      
        xOffs <- signif(diff(graphics::par("usr")[1:2])/170,3)
        yOffs <- signif(diff(graphics::par("usr")[3:4])/90,3)
        noNa <- if(is.null(names(tmP[useL2]))) 1:length(tmP) else which(is.na(names(tmP)[useL2]))
        if(length(noNa) >0 & all(annotColumn %in% colnames(merg))) names(tmP)[useL2[noNa]] <- merg[useLi[useL2[noNa]], wrMisc::naOmit(match(annotColumn[-1], colnames(merg)))[1]]
        if(length(NbestCol) <1) NbestCol <- 1
        displTx <- names(tmP[useL2])
        chNa <- is.na(displTx)
        if(any(chNa)) {displTx[which(chNa)] <- "unknown"; cexTxLab <- c(cexTxLab,cexTxLab*0.7)[1+chNa]}   # smaller label for 'unknown'
        if(all(chNa)) {if(!silent) message(fxNa," no names for display of best")
        } else graphics::text(Avalue[useLi[useL2]] +xOffs, Mvalue[useLi[useL2]] +yOffs, displTx, cex=cexTxLab, col=NbestCol, adj=0) }
    }              

    ## legend (if multiple symbols)
    pch[which(is.na(pch))] <- -2
    ch1 <- unique(pch)
    if(length(ch1) >1) {
      legInd <- which(!duplicated(merg[which(passAll), annotColumn[1]], fromLast=FALSE))
      legPch <- pch[which(passAll)[legInd]]
      legCol <- useCol[which(passAll)[legInd]]
      legBg <- ptBg[which(passAll)[legInd]]
      if(alph2 <1) {legCol <- substr(legCol,1,7); legBg <- substr(legBg,1,7)}  # reset to no transparency
      legLab <- merg[which(passAll)[legInd], annotColumn[1]]
      chNa <- is.na(legLab)
      if(any(chNa)) legLab[chNa] <- "NA"
      legOr <- if(length(legLab) >1) order(legLab) else 1   # not used so far 
      legLoc <- checkForLegLoc(merg[which(passAll),c("Avalue","Mvalue")] , sampleGrp=legLab, showLegend=FALSE)
      legCex <- stats::median(c(useCex,cexTxLab,1.2), na.rm=TRUE)
      graphics::legend(legLoc$loc, legend=legLab, col=if(FDR4color) grDevices::grey(0.5) else legCol, text.col=1, pch=legPch, if(length(ptBg) >0) pt.bg=ptBg, cex=legCex, pt.cex=1.2*legCex, xjust=0.5, yjust=0.5)  # as points
    }
  ## export results
  if(returnData) {
    merg <- merg[,-1*c(1,ncol(merg))]        # remove col 'ID' 'redundant' & 'pch'
    annCo <- wrMisc::naOmit(match(annotColumn, colnames(merg)))
    if(length(annCo) >0) cbind(merg[,annCo],  merg[,-annCo]) else merg }
  } }
     
