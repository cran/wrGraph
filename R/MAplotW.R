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
#' @param callFrom (character) allow easier tracking of messages produced
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging 
#' @return This function plots an MA-plot (to the current graphical device); if \code{returnData=TRUE}, a data.frame with ($ID, $Mvalue, $Avalue, $FDRvalue, $passFilt) gets returned  
#' @seealso (for PCA) \code{\link{plotPCAw}}
#' @examples
#' library(wrMisc)
#' set.seed(2005); mat <- matrix(round(runif(600),2), ncol=6)
#' rownames(mat) <- paste(rep(letters[1:25],each=4), letters[2:26])
#' MAplotW(mat[,2] -mat[,1], A=rowMeans(mat[,1:2]))
#' 
#' ## assume 2 groups with 3 samples each
#' matMeans <- rowGrpMeans(mat, gr=gl(2,3,labels=LETTERS[3:4]))
#' MAplotW(M=matMeans[,2] -matMeans[,1], A=rowMeans(mat)) 
#' 
#' ## assume 2 groups with 3 samples each and run moderated t-test (from package 'limma')
#' tRes <- moderTest2grp(mat, gl(2,3))
#' MAplotW(tRes$Mval, tRes$Amean)                          
#' MAplotW(M=tRes$Mval, A=rowMeans(tRes$means), FCth=1.3) 
#' MAplotW(tRes)
#' MAplotW(tRes, limM=c(-2,2), FCth=1.3) 
#' 
#' @export
MAplotW <- function(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, subTxt=NULL,
  grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limA=NULL,
  annotColumn=c("SpecType","GeneName","EntryName","Accession","Species","Contam"), annColor=NULL, cexPt=NULL, cexSub=NULL, cexTxLab=0.7, 
  namesNBest=NULL, NbestCol=1, NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE, debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="MAplotW")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  opar <- graphics::par(no.readonly=TRUE) 
  opar2 <- opar[-which(names(opar) %in% c("fig","fin","font","mfcol","mfg","mfrow","oma","omd","omi"))]    #
  on.exit(graphics::par(opar2))     # progression ok
  ## during function changes in  $mar,$cex.main,$cex.lab,$las 
  plotByFdr <- TRUE  #
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(Avalue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)           # grey
  fcRGB <- c(1,0,0)                  # red        for points passing  FC filt line
  multiComp <- TRUE                    # initialize  
  splNa <- annot <- ptType <- colPass <- ptBg <- grpMeans <- pcol <- pwComb <- FDRvalue <- FDRcol <- FdrList <- FDRty <- useComp <- useCompNa <- allCompNa <- pwSep <- pwNames <- NULL      # initialize
  if(length(pch) <1) pch <- 16
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  if(debug) message(fxNa,"Length Mvalue ",length(Mvalue)," ; length Avalue ",length(Avalue)," ; useComp ",useComp)
  if(identical(col,"FDR")) {FDR4color <- TRUE; col <- NULL} else FDR4color <- FALSE
  if(length(Mvalue) <1) message("Nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if(length(cexTxLab) <0) cexTxLab <- 0.7
    if("MArrayLM" %in% class(Mvalue) || "list" %in% class(Mvalue)) {
      if(debug) {message(fxNa,"maP0"); maP0 <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,fxNa=fxNa )}
      ## try working based on MArrayLM-object (Mvalue)
      if(debug) message(fxNa," '",namesIn[1],"' is list or MArrayLM-object ")
      ## initial check of useComp
      if(length(useComp) >1) { useComp <- wrMisc::naOmit(useComp)[1]
        if(!silent) message(fxNa,"Argument 'useComp' should be integer of length=1; using only 1st entry") }
      if(length(useComp) <1) { useComp <- 1
        if(!silent) message(fxNa,"Argument 'useComp' invalid, setting to 1") }
        
        
      ## address multiple questions (via useComp): need to check 'useComp', thus need to locate FDRvalues or pValues for group-names ..
      ## extract data: find suitable columns
      pcol <- wrMisc::naOmit(match(c("p.value","pvalue","pval","p"), tolower(names(Mvalue))))
      if(length(pcol) >0) names(pcol) <- c("p.value","pvalue","pval","p")[pcol]          
      FDRcol <- wrMisc::naOmit(match(c("fdr","bh","lfdr","by","bonferroni"), tolower(names(Mvalue))))
      if(length(FDRcol) >0) names(FDRcol) <- c("FDR","BH","lfdr","BY","bonferroni")[FDRcol]          
      if(debug) {message(fxNa,"maP1"); maP1 <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol )}
      
      if(length(pcol) ==0) stop(fxNa,"Cannot find suitable element for p-values")
      if(length(FDRcol) ==0) message(fxNa,"Cannot find suitable element for FDR-values")
      
      if("means" %in% names(Mvalue)) grpMeans <- Mvalue$means else stop(fxNa,"Need elment 'means' (group-means) in object 'Mvalue'")
      ## need to know which of pairwise comparisons get addressed:  need to understand setup -either from character useComp  or index & p-values
      if(debug) {message(fxNa,"maP1b"); maP1b <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans )}

      if("design" %in% names(Mvalue) ) {
        if(length(Mvalue$design)==0 || length(dim(Mvalue$design)) !=2 || nrow(Mvalue$design) ==0) stop(fxNa,"Unable to locate element with names of pairwise groups !")
        if(colnames(Mvalue$design)[1] =="(Intercept)" && ncol(Mvalue$design) ==2) {
          ## single comparison case (colnames of Mvalue$design are fixed); determine optimal sep and create pairwise name(s)
          grp <- colnames(Mvalue$means)            # 
          if(utils::packageVersion("wrMisc") < "2.0.0") {
            pwSep <- " "
          } else {
            pwSep <- " "
            #pwSep <- if(length(Mvalue$setup$sep)==1) Mvalue$setup$sep else wrMisc::getPWseparator(grp=grp)
          }
          pwNames <- if(length(Mvalue$setup$pwNames)==1) Mvalue$setup$pwNames else utils::combn(grp, 2)
          allCompNa <- if(length(Mvalue$setup$allCompNa) >0) Mvalue$setup$allCompNa else paste(pwNames[1,], pwNames[2,], sep=pwSep)
        } else {
          grp <- colnames(Mvalue$design)
          if(utils::packageVersion("wrMisc") < "2.0.0") {
            pwSep <- " "
          } else {  
            pwSep <- " "
            #pwSep <- if(length(Mvalue$setup$sep)==1) Mvalue$setup$sep else wrMisc::getPWseparator(grp=colnames(Mvalue$design), includeGrp=FALSE, silent=silent, debug=debug, callFrom=fxNa) 
          }
          allCompNa <- if(length(Mvalue$setup$allCompNa) >0) Mvalue$setup$allCompNa else colnames(Mvalue[[c(pcol,FDRcol)[1] ]])
        }

        ## investigate pairwise names from testing (& get separator) - need to find out which element of Mvalue
        if(length(allCompNa) ==0 ) stop(fxNa,"Unable to locate element with names of pairwise groups !")
        if(debug) {message(fxNa,"maP1c"); maP1c <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,grp=grp,pwSep=pwSep,pwNames=pwNames,allCompNa=allCompNa,useComp=useComp,useCompNa=useCompNa )}

        if(length(useCompNa) ==0) useCompNa <- allCompNa[useComp] 
        names(useComp) <- useCompNa
        
        useCompNaS <- if(length(useCompNa) ==0 && length(pwNames) >0 ) t(pwNames)[useComp,] else unlist(strsplit(useCompNa, pwSep))  
        avInd <- match(useCompNaS, colnames(grpMeans))
        Mvalue$Mval <- grpMeans[,avInd[1]] - grpMeans[,avInd[2]] 
        if(debug) {message(fxNa,"maP1c"); maP1c <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,pwSep=pwSep,useCompNaS=useCompNaS,avInd=avInd,annotColumn=annotColumn )}

        useCompNaS <- unlist(strsplit(useCompNa, pwSep))  
        avInd <- match(useCompNaS, colnames(grpMeans))
        Mvalue$Mval <- grpMeans[,avInd[1]] - grpMeans[,avInd[2]] 
        if(length(Avalue) ==0) {
          ## try to construct Avalue (and M-value) : 
          Avalue <- (grpMeans[,avInd[1]] + grpMeans[,avInd[2]]) /2
        }
        if(debug) {message(fxNa,"maP2"); maP2 <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep )}
      } else stop(fxNa,"Argument 'Mvalue' is missing key element $design")

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
          if(length(pch) < length(Avalue) && length(unique(wrMisc::naOmit(ptType))) >1) {
            if(length(pch) >1 && !silent) message(fxNa," (invalid pch) using default 'pch' oriented by $annot and starting from 15")
            pch <- 14 + as.integer(as.factor(ptType))            
          } 
          useAnnCol <- wrMisc::naOmit(useAnnCol) }
        annot <- Mvalue$annot[,useAnnCol] 
        if(annotColumn[1] %in% colnames(annot)) annot[,annotColumn[1]] <- ptType   # update with NAs trasformed to "NA"
        
        ## filtering
        if("filter" %in% names(Mvalue) && ((length(filtFin)==1 && isTRUE(filtFin)) || length(filtFin)==0)) {
          filtFin <- Mvalue$filter[,which(colnames(Mvalue$filter)==useCompNa)]           
        } else {    ## check custom filtFin
          chFi <- is.logical(filtFin) && length(filtFin)==length(Avalue)
          if(!chFi) { filtFin <- NULL
            if(!silent) message(fxNa,"Argument 'filtFin' seems invalud (should be logical vector of length of number of elements for plotting )")
          } else if(debug) message(fxNa,"Using custom filtFin")
        }
        if(length(filtFin) >0) {
          Avalue[which(filtFin)] <- NA    # better to set NA so that annotation, filtering etc stays same
          Mvalue$Mavl[which(filtFin)] <- NA                   
        }
      }           
      Mvalue <- Mvalue$Mval                             # DISMISS rest of object !!! .....
          
    } else {
      ## regular values
      if(is.data.frame(Mvalue)) { rowNaM <- rownames(Mvalue) 
        Mvalue <- as.numeric(as.matrix(Mvalue))
        if(is.null(names(Mvalue)) && !is.null(rowNaM)) names(Mvalue) <- rowNaM
      }   
      if(is.data.frame(Avalue)) { rowNaA <- rownames(Avalue) 
        Avalue <- as.numeric(as.matrix(Avalue))
        if(is.null(names(Avalue)) && !is.null(rowNaA)) names(Avalue) <- rowNaA
      }
      if(length(dim(Avalue))==2) if(ncol(Avalue) >1) { if(!silent) message(fxNa,"Note : 'Avalue' contains ",ncol(Avalue)," columns; but only a SINGLE column is expected ... taking mean of them")
        Avalue <- rowMeans(Avalue) } else Avalue <- try(as.numeric(if(is.data.frame(Avalue)) as.matrix(Avalue) else Avalue))
      if(debug) {message(fxNa,"maP2b"); maP2b <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep)}
      if(length(Mvalue) != length(Avalue)) stop(fxNa,"Length of Mvalue and Avalue do NOT match")
      if(length(names(Mvalue)) >0 && length(names(Avalue)) >0) {
        chNa <- names(Mvalue) != names(Avalue) 
        if(any(chNa)) warning(fxNa,"Note :  ",sum(chNa)," names (out of ",length(chNa),")  from Mvalue and Avalue do NOT match !!")
      }
      if(length(filtFin)==length(Mvalue) && is.logical(filtFin) && any(!filtFin)) {
        if(all(filtFin)) warning(fxNa,"All data get filtered away !!")
        Avalue[which(filtFin)] <- NA    # better to set NA so that annotation, filtering etc stays same
        Mvalue[which(filtFin)] <- NA        
      }
      if(debug) {message(fxNa,"maP2c"); maP2c <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep, pch=pch,annot=annot )}
      ## check order
      if(length(names(Mvalue)) >0 && length(names(Mvalue))==length(names(Mvalue))) {
        if(identical(names(Mvalue), names(Avalue))) {
          chOrd <- match(names(Avalue), names(Mvalue))
          if(any(is.na(chOrd))) warning(fxNa,"Caution : Unable to MATCH NAMES of Mvalue and Avalue; plot may be incorrect") else {
            if(!all(chOrd==1:length(Mvalue))) { Avalue <- Avalue[chOrd] 
              if(!silent) message(fxNa,"Note : ORDER of 'Mvalue' and 'Avalue' has been adjusted to 'Mvalue' (caution : 'filtFin' used in order of 'Mvalue' )") }
      } } }
      grpMeans <- rep(NA, length(Mvalue))    # could one provide from outside ? (useful ?)
    }
    if(debug) {message(fxNa,"maP3"); maP3 <- list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep, pch=pch,annot=annot )}

    ## prepare for plotting      
          
    if(length(pch)==1 && length(Mvalue) >1) pch <- rep(pch, length(Mvalue))
    if(length(pch) != length(Mvalue) && (length(pch) >1)) { if(!silent) message(fxNa,"Bizzare entry for 'pch'")
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
    if(debug) {message(fxNa,"maP8"); maP3b <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol )}

    ## adjust col (color) & pch
    if(!any(c(1,length(Mvalue)) %in% length(pch))) {
      if(!silent) message(fxNa,"Argument 'pch' should be either length=1 or correspond to length of data, reset to default=16")
      pch <- 16 }
    if(length(col) >1 && length(col) <length(Mvalue)) {
      if(!silent) message(fxNa,"Argument 'col' should be either length=1 or correspond to length of data, reset to default=NULL")
      col <- NULL }

    if(debug) message(fxNa," ++ DONE extracting columns : ",wrMisc::pasteC(colnames(merg),quo="'"))

    ## apply filtering (keep all lines where at least one condition passes)
    if(!identical(filtFin, FALSE) && length(filtFin)==nrow(merg)) {                         #  use filtering provided
      if(sum(filtFin) >0 && sum(filtFin) < nrow(merg)) { 
        whFilt <- which(merg$filtFin)
        if(length(pch) >1) pch <- pch[whFilt]
        if(length(col) >1) col <- col[whFilt]
        merg <- merg[whFilt,]
        if(!silent) message(fxNa,"Filtered (based on 'filtFin') from ",length(filtFin)," to  ",nrow(merg)," lines")
      }
    } else filtFin <- rep(TRUE, nrow(merg))
    if(debug) {message(fxNa,"maP4"); maP4 <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep )}
    
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
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","MA-Plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "MA-plot"
    if(length(FCthrs) <1) FCthrs <- 1.5 
    #needed?#if(length(FdrThrs) <1) FdrThrs <- 0.05 
    if(debug) {message(fxNa,"maP4b"); maP4b <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep )}

    ## count no of passing
    passAll <- passFC <- if(length(FCthrs) ==1 && !any(is.na(FCthrs))) abs(merg[,"Mvalue"]) >= log2(FCthrs) else merg[,"filtFin"]      ## convert FCthrs to log2
    passAll <- merg[,"filtFin"] & passFC    

    chNA <- is.na(passAll)                              # passFdr may contain NAs
    if(any(chNA)) passAll[which(chNA)] <- FALSE 
    if(debug) message(fxNa,"  ",sum(passFC,na.rm=TRUE)," passing FCthrs ")
    ## color for points passing filter
    if(length(col) >0) if(length(col) != nrow(merg)) { col <- NULL
      if(!silent) message(fxNa,"Invalid entry for 'col', should be of length=",nrow(merg),", resetting to default")}
    if(debug) {message(fxNa,"maP5"); maP5 <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep,passAll=passAll,col=col,grayIncrem=grayIncrem,cexPt=cexPt,annotColumn=annotColumn,annColor=annColor )}

    if(length(col) <1) {
      alph <- sort(c(0.14, round(0.6/log10(length(Mvalue)),2), 0.8))[2]       # alph <- round(12/sqrt(nrow(eBayesLst$pValue)),2)
      alph2 <- sort(c(round(7/(5 +sum(passAll)^0.7),2), alph,0.9))[2]                   # for points passing thresholds
      useCol <- if(grayIncrem) grDevices::rgb(0.35,0.35,0.35,alph) else grDevices::rgb(0.7,0.7,0.7)  # basic color (grey70)
      useCex <- if(length(cexPt) >0) cexPt else max(round(0.8 + 2/(1 +sum(filtFin, na.rm=TRUE))^0.28,2), 1.1)
      
      chCol <- unique(merg[, annotColumn[1]])      # check how many different colors may be needed
      chNaC <- is.na(chCol)
      if(any(chNaC)) chCol[which(chNaC)] <- "NA"          # transform NA to 'NA'
      if(length(annColor) >0) {colPass <- annColor} else { if(length(chCol) >4) {
        colPass <- cbind(red=c(141,72,90,171, 220,253,244,255), green=c(129,153,194,221, 216,174,109,0), blue=c(194,203,185,164, 83,97,67,0))       
        colPass <- grDevices::rgb(red=colPass[,1], green=colPass[,2], blue=colPass[,3], alph2, maxColorValue=255)
        if(length(chCol) >8) { colPass <- c(colPass, rep(colPass[8], length(chCol) -8))
          if(!silent) message(fxNa," > 8 different groups found, using 8th color after 7th group")}
      } else colPass <- grDevices::rgb(c(0.95,0.2,0,0.75), c(0.15,0.2,0.9,0.35), c(0.15,0.95,0,0.8), alph2)}    # red, blue, green, purple (luminosity adjusted) 
      useCol <- rep(useCol[1], nrow(merg))         # fuse basic gray to colors for different types

      ## integrate names of annColor as order of colPass
      if(length(names(annColor)) >0) {
        uniTy <- unique(merg[which(passAll),annotColumn[1]])
        colPass <- colPass[match(names(annColor), uniTy)]
      }     
      ## assign color for those passing
      if(any(passAll)) { if(FDR4color) {
          useCol[which(passAll)] <- .colorByPvalue(merg$FDR[which(passAll)])          # ultimately replace/integrate to wrMisc::colorAccording2()
        } else {
        useCol[which(passAll)] <- colPass[ if(length(unique(merg[which(passAll),annotColumn[1]])) >1) wrMisc::levIndex(
          merg[which(passAll),annotColumn[1]]) else rep(1,sum(passAll))] } }  # assign colors for those passing 
    } else useCol <- col
    if(debug) {message(fxNa,"maP6"); maP6 <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,useCompNa=useCompNa,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep,namesNBest=namesNBest,passAll=passAll,useCol=useCol )}
        
    ## adjust fill color for open symbols
    chPch <- pch %in% c(21:25)
    if(any(chPch)) { ptBg <- useCol
      ptBg[which(chPch)] <- useCol[which(chPch)]    # background color for filled symbols
      useCol[which(chPch)] <- 1                     # contour as black
    }

    ## main graphic
    tmp <- try(graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa, las=1), silent=TRUE)
    ## rather directly plot FDR
    tmp <- try(graphics::plot(merg[,"Avalue"], merg[,"Mvalue"], pch=pch, cex=useCex, main=tit1, 
      ylab="M value (log FC)", col=useCol, xlab=xLab, cex.lab=cexLa, xlim=limM,ylim=limA, pt.bg=ptBg), silent=TRUE)
    if(inherits(tmp, "try-error")) warning(fxNa,"UNABLE to produce plot !") else {
      sTxt <- if(length(subTxt) ==1) subTxt else { if(multiComp) paste0(if(length(names(useComp)) >0) names(useComp) else NULL)}
      #old# sTxt <- if(length(subTxt) ==1) subTxt else { if(multiComp) paste0(if(length(names(useComp)) >0) names(useComp) else paste0("useComp=",useCompNa," (",useComp,")"),"; ",collapse="")}
      sTxt <- paste0(sTxt," n=",length(Mvalue),
        if(!all(is.na(c(FCthrs)))) paste("; ",sum(passAll, na.rm=TRUE),"(color) points passing", if(!is.na(FCthrs)) paste0("FCthr=", as.character(signif(FCthrs,3))) ))
      graphics::mtext(sTxt, cex=0.75, line=0.2)    
      
      if(!all(is.na(FCthrs))) { 
        if(debug) message(fxNa," n=",length(Mvalue),"  FCthrs=",as.character(FCthrs),"  filt.ini=", sum(filtFin, na.rm=TRUE),"  passAll=",sum(passAll,na.rm=TRUE),
          " ; range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=TRUE),3))," ;  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
        graphics::abline(h=c(-1,1)*(log2(FCthrs) + diff(graphics::par("usr")[1:2])/500), col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
      
      if(debug) {message(fxNa,"maP7"); maP7 <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep,namesNBest=namesNBest,passAll=passAll )}
      
      ## add names to best points
      if(length(namesNBest) >0) { 
        if(any(sapply( c("passThr","pass","passFC"), identical, namesNBest))) namesNBest <- sum(passAll)
        if(!is.integer(namesNBest)) namesNBest <- try(as.integer(namesNBest), silent=TRUE)
        if(namesNBest >0 && any(passAll)) {      
          useLi <- if(any(!passAll)) which(passAll) else 1:nrow(merg)
          tmP <- as.numeric(merg[useLi,"Mvalue"])
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
          if(length(noNa) >0 && all(annotColumn %in% colnames(merg))) names(tmP)[useL2[noNa]] <- merg[useLi[useL2[noNa]], wrMisc::naOmit(match(annotColumn[-1], colnames(merg)))[1]]
          if(length(NbestCol) <1) NbestCol <- 1
          displTx <- names(tmP[useL2])
          chNa <- is.na(displTx)
          if(any(chNa)) {displTx[which(chNa)] <- "unknown"; cexTxLab <- c(cexTxLab,cexTxLab*0.7)[1+chNa]}   # smaller label for 'unknown'
          if(debug) {message(fxNa,"maP7b"); maP7b <- list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,FDRcol=FDRcol,grpMeans=grpMeans,allCompNa=allCompNa,pwSep=pwSep,namesNBest=namesNBest,passAll=passAll,chNa=chNa,useLi=useLi,useL2=useL2,xOffs=xOffs,yOffs=yOffs,cexTxLab=cexTxLab,NbestCol=NbestCol)}
          if(all(chNa)) {if(!silent) message(fxNa,"No names for display of best")
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
      } }
  tmp <- try(graphics::par(mar=opar$mar, cex.main=opar$cex.main, las=opar$las), silent=TRUE)
      
  ## export results
  if(returnData) {
    merg <- merg[,-1*c(1,ncol(merg))]        # remove col 'ID' 'redundant' & 'pch'
    annCo <- wrMisc::naOmit(match(annotColumn, colnames(merg)))
    if(length(annCo) >0) cbind(merg[,annCo],  merg[,-annCo]) else merg }
  } }
        
