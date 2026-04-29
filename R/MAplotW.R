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
#' @param Avalue (numeric, list or data.frame) if \code{NULL} it is assumed that A-values can be extracted form argument \code{Mvalue}
#' @param useComp (integer) choice of one of multiple comparisons present in \code{Mvalue} (if generated using \code{moderTestXgrp()})  
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied
#' @param tit (character) custom title
#' @param ProjNa (character) add project-name to (automatic) title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC)
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points  (not exceeding thresholds) as increased shades of gray
#' @param col (character) custom color(s) for points of plot (see also \code{\link[graphics]{par}})
#' @param pch (integer) type of symbol(s) to plot (default=16) (see also \code{\link[graphics]{par}})
#' @param compNa (character) names of groups compared
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace
#' @param cexMa (numeric) font-size of title, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexLa (numeric) size of axis-labels, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param limM (numeric, length=2) range of axis M-values
#' @param limp (numeric, length=2) range of axis FDR / p-values
#' @param annotColumn (character) column names of annotation to be extracted (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot).
#'   The first entry (typically 'SpecType') is used for different symbols in figure, the second (typically 'GeneName') is used as prefered text for annotating the best points (if \code{namesNBest} allows to do so.)
#' @param annColor (character or integer) colors for specific groups of annotation (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param expFCarrow (logical, character or numeric) optional adding arrow for expected fold-change; if \code{TRUE} the expected ratio will be extracted from numeric concentration-indications from sample-names
#'  if \code{numeric} an arrow will be drawn (M-value as 1st position, color of 2nd position of vector).
#' @param cexPt (numeric) size of points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexSub (numeric) size of subtitle, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexTxLab (numeric) size of text-labels for points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param namesNBest (integer or character) for display of labels to points in figure: if 'pass','passThr' or 'signif' all points passing thresholds; if numeric (length=1) this number of best points will get labels
#'   if the initial object \code{Mvalue} contains a list-element called 'annot' the second of the column specified in argument \code{annotColumn} will be used as text
#' @param NbestCol (character or integer) colors for text-labels of best points, also used for arrow
#' @param colBySpecType (logical) incase arument \code{Mvalue} is MArrayLM-object it is possible to use different color-codes for points passing thresholds based on Mvalue$annot[,"SpecType"]; use this with multi-species benchmark tests 
#' @param sortLeg (character) sorting of 'SpecType' annotation either ascending ('ascend') or descending ('descend'), no sorting if \code{NULL}
#' @param NaSpecTypeAsContam (logical) consider lines/proteins with \code{NA} in Mvalue$annot[,"SpecType"] as contaminants (if a 'SpecType' for contaminants already exits)
#' @param useMar (numeric, length=4) custom margings (see also \code{\link[graphics]{par}})
#' @param returnData (logical) optional returning data.frame with (ID, Mvalue, pValue, FDRvalue, passFilt)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @param debug (logical) additional messages for debugging
#' @return This function plots an MA-plot (to the current graphical device); if \code{returnData=TRUE}, a data.frame with ($ID, $Mvalue, $Avalue, $passFilt) gets returned  
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
MAplotW <- function(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, tit=NULL, ProjNa=NULL, FCthrs=NULL,
  subTxt=NULL, grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limp=NULL,
  annotColumn=c("SpecType","GeneName","EntryName","Accession","Species","Contam"), annColor=NULL, expFCarrow=FALSE,cexPt=NULL, cexSub=NULL,
  cexTxLab=0.7, namesNBest=NULL, NbestCol=1, colBySpecType=FALSE, sortLeg="descend", NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE, debug=FALSE) {

#MAplotW <- function(Mvalue, Avalue=NULL, useComp=1, filtFin=NULL, tit=NULL, ProjNa=NULL, FCthrs=NULL, subTxt=NULL,
#  grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limA=NULL,
#  annotColumn=c("SpecType","GeneName","EntryName","Accession","Species","Contam"), annColor=NULL, cexPt=NULL, cexSub=NULL, cexTxLab=0.7, 
#  namesNBest=NULL, NbestCol=1, NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE, debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="MAplotW")
  opar <- graphics::par(no.readonly=TRUE)
  opar2 <- opar[-which(names(opar) %in% c("fig","fin","font","mfcol","mfg","mfrow","oma","omd","omi"))]    #
  on.exit(graphics::par(opar2))     # progression ok
  ## during function changes in  $mar,$cex.main,$cex.lab,$las 
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(Avalue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)           # grey
  fcRGB <- c(1,0,0)                  # red        for points passing  FC filt line
  inversedComp <- singleCompSetup <- FALSE                # inversedComp: if useComp is given as text it may be inverse of p-values given; singleCompSetup : cases sof single comparison BUT 2 cols of p.values
  splNa <- annot <- ptType <- colPass <- ptBg <- grpMeans <- pcol <- pwComb <- pwNames <- useCompNa <- sortLeg <- NULL      # initialize    # FDRvalue <- FdrList <- FDRty <-

  .extrSep2 <- function(comb , indiv) {     ## extract separator by striping indiv elements (iniv) from combined pairwise names (comb)  .. move to wrMisc ??
    fx2 <- function(x, y) sub(paste0(wrMisc::protectSpecChar(x),"$"),"", sub(paste0("^",wrMisc::protectSpecChar(x)),"", y))   # remove x from head and tail   (eg rm C from end)
    for(i in indiv) comb <- fx2(i, comb)
    unique(comb[which(nchar(comb) >0)])   # will/might return NULL when sep==""
  }

  if(length(pch) <1) pch <- 16
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  if(debug) message(fxNa,"Length Mvalue ",length(Mvalue)," ; length Avalue ",length(Avalue)," ; useComp ",useComp)
  #if(identical(col,"FDR")) {FDR4color <- TRUE; col <- NULL} else FDR4color <- FALSE
  if(length(Mvalue) <1) message("Nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if(utils::packageVersion("wrMisc") <= "2.0.2" && "MArrayLM" %in% class(Mvalue) || "list" %in% class(Mvalue))
    message("Please install latest version of wrMisc (> 2.0.1) from CRAN ")
    Mvalue <- NULL
  }  
  if(length(Mvalue) >0) {                   # main
    ## data seem valid to make MAplot
    if(length(cexTxLab) <0) cexTxLab <- 0.7
    if("MArrayLM" %in% class(Mvalue) || "list" %in% class(Mvalue)) {
      if(debug) { message(fxNa,"MAP0"); MAP0 <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,singleCompSetup=singleCompSetup,fxNa=fxNa )}
      ## extract setup from MArrayLM
      setup <- wrMisc::getPairwiseSetup(Mvalue, silent=silent, debug=debug, callFrom=fxNa)

      ## initial check of useComp & map (point to  alphabet sorted setup$grpNa)
      if(length(useComp) >1) useComp <- wrMisc::naOmit(useComp)
      if(length(useComp) >2) { useComp <- useComp[1:2]; warning(fxNa,"Argument 'useComp' seems invalid, trimming to 2 (2 groups to be compared)")}
      if(length(useComp) <1) { useComp <- 1
        if(!silent) message(fxNa,"Argument 'useComp' seems absent, setting to default =1") }
      if(is.character(useComp)) { 
        useComp2 <- if(length(useComp) ==2) which(setup$pwGrpNa[,1, drop=FALSE]==useComp[1] & setup$pwGrpNa[,2, drop=FALSE]==useComp[2]) else which(rownames(setup$index)==useComp) 
        if(length(useComp2) ==0) {
          ## could be inverse question, thus need the to reverse M-values on plot
          useComp2 <- if(length(useComp) ==2) which(setup$pwGrpNa[,1, drop=FALSE]==useComp[2] & setup$pwGrpNa[,2, drop=FALSE]==useComp[1]) else which(paste0(setup$pwGrpNa[,2],setup$sep,setup$pwGrpNa[,1])==useComp)   
          if(length(useComp2)==0) stop(fxNa,"UNABLE TO MAP character 'useComp' to experimental setup") else inversedComp <- TRUE
        } 
        useComp <- useComp2
      }
      if(debug) { message(fxNa,"MAP0b"); MAP0b <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,setup=setup,fxNa=fxNa )}
      if(is.numeric(useComp)) {
        if(any(useComp <1 | useComp > nrow(setup$pwGrpNa))) stop(fxNa,"Invalid entry for 'useComp' (if nuermic should be index of already tested groups)")  
        useComp <- as.integer(useComp) 
        if(length(useComp)==2) {
          useComp2 <- which(setup$index[,1, drop=FALSE]==useComp[1] & setup$index[,2, drop=FALSE]==useComp[2])
          if(length(useComp2) ==0) { useComp2 <- which(setup$index[,1, drop=FALSE]==useComp[2] & setup$index[,2, drop=FALSE]==useComp[1])    ## could be inverse question, thus need the to reverse M-values on plot
            if(length(useComp2) ==0) stop(fxNa,"Unable to MAP numeric 'useComp' to experimental setup") else inversedComp <- TRUE}
          useComp <- useComp2
        }
      }
      useCoPw <- setup$index[useComp, ]    # index as length=2 (for $means etc)
      if(inversedComp) { useCoPw <- useCoPw[2:1] 
        if(!silent) message(fxNa,"Note : Question is inversed to precalculated test (M-values have been adjusted)") }
      
      if(debug) {message(fxNa,"MAP1"); MAP1 <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,useCoPw=useCoPw,setup=setup )}
      
      if("means" %in% names(Mvalue)) grpMeans <- Mvalue$means else stop(fxNa,"Need element 'means' (group-means) in object 'Mvalue'")
      ## need to know which of pairwise comparisons get addressed:  need to understand setup -either from character useComp  or index & p-values
      if(debug) {message(fxNa,"MAP1b"); MAP1b <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,grpMeans=grpMeans,setup=setup,useCoPw=useCoPw )}

      
      Mvalue$Mval <- grpMeans[,useCoPw[1]] - grpMeans[,useCoPw[2]]
      ##++##
      Avalue <- rowMeans(grpMeans[,useCoPw[1:2]], na.rm=TRUE)
      #pValue <- if(length(pcol)==1) Mvalue[[pcol]][,useComp] else NULL
      if(debug) {message(fxNa,"MAP2"); MAP2 <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,grpMeans=grpMeans,annotColumn=annotColumn,useCoPw=useCoPw,setup=setup,colBySpecType=colBySpecType )}


      ## recuperate $annot if present and use for symbol
      if("annot" %in% names(Mvalue) && isTRUE(colBySpecType)) {
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
        if(debug) {message(fxNa,"MAP2b"); MAP2b <-  list(Mvalue=Mvalue,Avalue=Avalue,annot=annot,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,grpMeans=grpMeans,annotColumn=annotColumn,useAnnCol=useAnnCol,useCompNa=useCompNa,pch=pch,useCoPw=useCoPw )}   
      }  
        
      ## filtering
      if(("filter" %in% names(Mvalue)) && isTRUE(filtFin)) {
        # extract existing filtering
        filtFin <- Mvalue$filter[,useComp]           
        #filtFin <- Mvaue$filter[,which(colnames(Mvalue$filter)==useCompNa)]           
      } else {    ## check if custom filtFin
        chFi <- length(filtFin)==length(Avalue)
        if(!chFi) { 
          if(!isFALSE(filtFin) && !silent) message(fxNa,"Argument 'filtFin' seems invalid (should be TRUE or logical vector of length of number of elements for plotting )")
          filtFin <- NULL          
        } else if(debug) message(fxNa,"Custom filtFin seems valid")
      }
      if(length(filtFin) >0) {
        if(debug) message(fxNa,"Filtering : Setting ",sum(!filtFin & !is.na(Avalue))," additional values as NA")                   
        Avalue[which(!filtFin)] <- NA    # better to set NA so that annotation, filtering etc stays same 
        ## Mvalue : if only 2 groups $Mval should be OK
        if("Mval" %in% names(Mvalue)) Mvalue$Mval[which(!filtFin)] <- NA
      }
      if(debug) {message(fxNa,"MAP2c"); MAP2c <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,filtFin=filtFin,ProjNa=ProjNa,FCthrs=FCthrs,pcol=pcol,grpMeans=grpMeans,annotColumn=annotColumn,setup=setup )}
      Mvalue <- Mvalue$Mval                             # DISMISS rest of object !!! .....
      singleCompSetup <- !(length(dim(setup$pwGrpNa))==2 && nrow(setup$pwGrpNa) >1)
                
    } else {
      ## ... case of explicit Avalue argument
      singleCompSetup <- TRUE
      if(length(Avalue) <1) stop("Argument 'Avalue' is required (if 'Mvalue' not 'MArrayLM'-type object) !")
      if(length(dim(Avalue)) >1) if(ncol(Avalue) >1) {
        if(!silent) message(fxNa,"Note, ",namesIn[2]," has ",ncol(Avalue)," columns, using last column")
        pNa <- rownames(Avalue)
        Avalue <- as.numeric(Avalue[,ncol(Avalue)] )     # use last col if multiple
        names(Avalue) <- pNa }
      if(length(Avalue) != length(Mvalue)) warning(fxNa,"Length of Mvalue and Avalue DIFFER !!")  
      ## filtering
      chFi <- length(filtFin)==length(Avalue)
      if(!chFi) { 
        if(length(filtFin) >1 && !silent) message(fxNa,"Argument 'filtFin' seems invalid (should be logical vector of length of number of elements for plotting )")
        filtFin <- NULL          
      } else { if(debug) message(fxNa,"Custom filtFin seems valid")
        if(any(!filtFin)) Avalue[which(!filtFin)] <- NA    # better to set NA so that annotation, filtering etc stays same
      }
    }    
    if(debug) {message(fxNa,"MAP8"); MAP8 <-  list(Mvalue=Mvalue,Avalue=Avalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans,singleCompSetup=singleCompSetup)}


    ## need to introduce -log10 to Avalue
    chNA <- is.na(Avalue)
    if(all(chNA)) stop(fxNa,"All p-values are NA, nothing to draw !")
    if(debug && any(Avalue <0, na.rm=TRUE)) message(fxNa,"Some p-values are negative, this should not be !  Maybe log values were given by error ?")

    ## check for (same) order, adjust Mvalue & Avalue according to names
    chNa <- list(MNa=if(length(dim(Mvalue)) >1) rownames(Mvalue) else names(Mvalue),
      pNa=if(length(dim(Avalue)) >1) rownames(Avalue) else names(Avalue))
    nIni <- c(M=length(Mvalue), p=length(Avalue))
    if(length(chNa$MNa) >0 && length(chNa$pNa) >0) {        # ie both have names, so one can match names
      if(!all(chNa$MNa==chNa$pNa, na.rm=TRUE)) {
        matchNa <- wrMisc::naOmit(match(chNa$MNa, chNa$pNa))
        if(length(matchNa) <1) stop("Both 'Mvalue' and 'Avalue' have names, but none of them match !!")
        if(!all(matchNa, 1:length(Avalue), na.rm=TRUE)) {
          Avalue <- Avalue[matchNa]
          Mvalue <- wrMisc::naOmit(Mvalue[match(names(Avalue), names(Mvalue))]) }
      }
    } else {
      if(length(Mvalue) != length(Avalue)) stop("p- and M- values have different length, but no names to match !!  (M=",length(Mvalue)," vs p=",length(Avalue),")")
    }
    if(length(grpMeans) <1) grpMeans <- matrix(rep(NA,2*length(Mvalue)), ncol=2, dimnames=list(names(Mvalue),c("mean1","mean2")))
    if(debug) {message(fxNa,"MAP9"); MAP9 <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans, filtFin=filtFin,annotColumn=annotColumn,annot=annot, pch=pch,sortLeg=sortLeg,batchFig=batchFig,singleCompSetup=singleCompSetup)}       # pairwCol=pairwCol,

    ## prepare/integrate FILTERING
    if(length(filtFin) >0) {
      ## if filtFin is matrix use each line with min 1 instance of TRUE,
      if(length(dim(filtFin)) >1) filtFin <- as.logical(as.matrix(filtFin)[,useComp])    # use rows with >= 1 TRUE
      if(length(names(filtFin)) >0) {
        matchNa <- wrMisc::naOmit(match(rownames(Avalue), names(filtFin)))
        if(length(matchNa)==length(Avalue)) filtFin <- as.logical(filtFin[matchNa])
      }
    } else filtFin <- rep(TRUE, nrow(grpMeans))
    if(debug) {message(fxNa,"MAP9b"); MAP9b <-  list(Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans, filtFin=filtFin,annotColumn=annotColumn,annot=annot, pch=pch,sortLeg=sortLeg,batchFig=batchFig)}       # pairwCol=pairwCol,

    ## start creating merged data for plot (& export)
    merg <- if(length(annot) >0) data.frame(ID=NA, grpMeans, Mvalue=Mvalue, Avalue=as.numeric(Avalue),
      filtFin=filtFin, annot, pch=pch, stringsAsFactors=FALSE) else {
        data.frame(ID=NA, grpMeans, Mvalue=Mvalue, Avalue=as.numeric(Avalue), 
        filtFin=filtFin, pch=pch, stringsAsFactors=FALSE) }
    if(!"Mvalue" %in% colnames(merg)) colnames(merg)[match("Avalue", colnames(merg)) -1] <- "Mvalue"    
    if(length(names(Mvalue)) >0) merg[,1] <- names(Mvalue) else {if(length(names(Avalue)) >0) merg[,1] <- names(Avalue)}
    ## replace NA in 'SpecType' by 'NA'

    if(annotColumn[1] %in% colnames(merg)) { chNa <- is.na(merg[,annotColumn[1]])     # replace NAs in col "SpecType" by "NA"
      if(any(chNa)) merg[which(chNa), annotColumn[1]] <- "NA"
    } else { merg <- cbind(merg, rep(1,nrow(merg)))            # add colum for 'SpecType'
      colnames(merg)[ncol(merg)] <- annotColumn[1] }
    if(debug) {message(fxNa,"MAP10"); MAP10 <-  list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans, filtFin=filtFin,annotColumn=annotColumn,annot=annot, pch=pch,sortLeg=sortLeg,batchFig=batchFig)}
    #check#plot(Avalue ~Mvalue, data=MAP10$merg)   # ok

    ## adjust col & pch
    if(!any(c(1,length(Mvalue)) %in% length(pch))) {
      if(!silent) message(fxNa,"Argument 'pch' should be either length=1 or correspond to length of data, reset to default=16")
      pch <- 16 }
    if(length(col) >1 && length(col) <length(Mvalue)) {
      if(!silent) message(fxNa,"Argument 'col' should be either length=1 or correspond to length of data, reset to default=NULL")
      col <- NULL }
    if(debug) message(fxNa," ++ DONE extracting columns : ",wrMisc::pasteC(colnames(merg),quo="'"))

    ## apply filtering
    msg <- "Data provided in 'Mvalue' and 'Avalue' "
    if(debug) {message(fxNa,"MAP10")}

    if(!silent && nrow(merg) < round(length(Mvalue)/10)) message(fxNa," .. note : less than 10% of",msg," were matched") else {
      if(!silent && nrow(merg) < length(Mvalue)/2) message(fxNa," .. NOTE : less than 50% of",msg," were matched !!")}
    if(debug) message(fxNa,msg," were matched to ",nrow(merg)," common entries")
    ## apply filtering (keep all lines where at least one condition passes)
    if(length(filtFin) >0 && any(isFALSE(filtFin), na.rm=TRUE)) {                         #  use filtering provided
      if(sum(filtFin) >0 && sum(filtFin) < nrow(merg)) {
        whFilt <- which(merg$filtFin)
        if(length(pch) >1) pch <- pch[whFilt]
        if(length(col) >1) col <- col[whFilt]
        merg <- merg[whFilt,]
        if(!silent) message(fxNa,"Filtered (based on 'filtFin') from ",length(filtFin)," to  ",nrow(merg)," lines")
      }
    } else filtFin <- rep(TRUE, nrow(merg))
    if(debug) {message(fxNa,"MAP11"); MAP11 <-  list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans, filtFin=filtFin,annotColumn=annotColumn,annot=annot, pch=pch,sortLeg=sortLeg,batchFig=batchFig)}
    #check#plot(Avalue ~Mvalue, data=MAP11$merg)   # ok

    ## sort merg, so that legend always gets constructed the same order, ascending ('ascend') or descending ('descend')
    sortLeg <- if(identical(sortLeg,"ascend")) FALSE else {if(identical(sortLeg,"descend")) TRUE else NULL}
    if(length(sortLeg) ==1 && annotColumn[1] %in% colnames(merg)) merg <- merg[order(merg[,annotColumn[1]], decreasing=sortLeg),]

    ## update ..
    nIDco <- sum(c("ID","nredID","uniqID") %in% colnames(merg))                   #  number of heading columns in 'merg'
    Mvalue <- as.numeric(if("Mvalue" %in% colnames(merg)) merg[,"Mvalue"] else merg[,nIDco+1])
    Avalue <- as.numeric(if("Avalue" %in% colnames(merg)) merg[,"Avalue"] else {
      if(length(dim(Mvalue)) >0) merg[,ncol(Mvalue) +nIDco +1] else merg[,nIDco+2]})
    pch <- merg[,"pch"]            # update
    ptType <- if(annotColumn[1] %in% colnames(merg)) merg[,annotColumn[1]] else rep(1,nrow(merg))     # update "SpecType"
    if(debug) {message(fxNa,"MAP12"); MAP12 <- list(merg=merg,Mvalue=Mvalue,filtFin=filtFin)}

    ## prepare for  plotting
    if(is.null(cexSub)) cexSub <- cexLa +0.05
    yLab <- "M-value (log2 fold-change)"
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","MA-Plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "MA-plot"

    tit1 <- if(length(tit)==1 && nchar(tit) >0) tit else paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","MA-Plot"),
      if(!is.null(compNa)) {if(length(compNa)==1) compNa else c(compNa[1]," vs ",compNa[2])}), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "MA-plot"
    if(debug) {message(fxNa,"MAP12b"); MAP12b <-  list(merg=merg,Mvalue=Mvalue,Avalue=Avalue,useComp=useComp,grpMeans=grpMeans, filtFin=filtFin,annotColumn=annotColumn,annot=annot, pch=pch,sortLeg=sortLeg,batchFig=batchFig,FCthrs=FCthrs)}

    ## count no of passing
    if(length(FCthrs) != 0) {
      msg1 <- "Invalid entry for argument 'FCthrs' !  Must be single numeric and >1 (as lin scale fold-change limit)  .. ignoring" 
      chNum <- try(as.numeric(FCthrs))
      if(inherits(chNum, "try-error")) { FCthrs <- NULL
        if(!silent) message(fxNa,msg1)
      } else if(any(is.na(FCthrs)) || length(FCthrs) >1) FCthrs <- wrMisc::naOmit(FCthrs)[1]
      if(length(FCthrs) ==1 && FCthrs <= 1) { FCthrs <- NULL
        if(!silent) message(fxNa,msg1) }
    } 

    passFC <- if(length(FCthrs) ==1 && !any(is.na(FCthrs))) abs(merg[,"Mvalue"]) >= log2(FCthrs) else merg[,"filtFin"]      ## convert FCthrs to log2
    passAll <- merg[,"filtFin"] & passFC
    if(debug) {message(fxNa,"MAP12c")}

    chNA <- is.na(passAll)                              # 
    if(any(chNA)) passAll[which(chNA)] <- FALSE
    if(debug) message(fxNa,"  ",sum(passFC,na.rm=TRUE)," passing FCthrs")
    if(debug) {message(fxNa,"MAP13"); MAP13 <-  list(merg=merg,Mvalue=Mvalue,filtFin=filtFin,col=col,passAll=passAll,passFC=passFC,FCthrs=FCthrs,cexPt=cexPt,annotColumn=annotColumn,annColor=annColor)}
    #check#plot(Avalue ~Mvalue, data=MAP13$merg)   # ok

    ## color for points passing filter
    if(length(col) >0) if(length(col) != nrow(merg)) { col <- NULL
      if(!silent) message(fxNa," invalid entry for 'col', should be of length=",nrow(merg),", resetting to default")}
    if(length(col) <1) {
      alph <- sort(c(0.14, round(0.6/log10(length(Mvalue)),2), 0.8))[2]       # alph <- round(12/sqrt(nrow(eBayesLst$Avalue)),2)
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
        uniTy <- unique(merg[which(passAll), annotColumn[1]])
        colPass <- colPass[match(names(annColor), uniTy)]
      }
      if(debug) {message(fxNa,"MAP13a"); MAP13a <-  list(useCol=useCol,colPass=colPass, merg=merg,Mvalue=Mvalue,filtFin=filtFin,col=col,passAll=passAll,cexPt=cexPt,annotColumn=annotColumn,annColor=annColor)}

      ## assign color for those passing (only if FCthrs was defined)
      if(length(FCthrs) ==1) useCol[which(passAll)] <- if(any(passAll, na.rm=TRUE)) colPass[if(length(unique(merg[which(passAll),annotColumn[1]])) >1) wrMisc::levIndex(merg[which(passAll),annotColumn[1]]) else rep(1,sum(passAll))]  # assign colors for those passing
    } else useCol <- col
    if(debug) {message(fxNa,"MAP13b"); MAP13b <-  list(merg=merg,Mvalue=Mvalue,filtFin=filtFin,col=col,passAll=passAll,cexPt=cexPt,annotColumn=annotColumn,annColor=annColor,pch=pch,useCol=useCol)}
#}}  #ok

    ## adjust fill color for open symbols
    chPch <- pch %in% c(21:25)
    if(any(chPch, na.rm=TRUE)) { ptBg <- useCol
      ptBg[which(chPch)] <- useCol[which(chPch)]    # background color for filled symbols
      useCol[which(chPch)] <- 1                     # contour as black
    }
    ## define limA (equiv to limp)
    ## +++ 
    limA <- range(merg$Avalue, na.rm=TRUE)


  ## main graphics
  if(length(Mvalue) >0) {
    pl1 <- try(graphics::par(mar=if(length(useMar)==4) useMar else c(6.5,4,4,2), cex.main=cexMa, las=1), silent=TRUE)
    if(inherits(pl1, "try-error")) {Mvalue <- NULL; message("UNABLE TO SET PLOT MARGINS !!  check plotting device ...")}
    if(debug) {message(fxNa,"MAP14a"); MAP14a <-  list(pl1=pl1,merg=merg,Mvalue=Mvalue,filtFin=filtFin,useMar=useMar,FCthrs=FCthrs,pch=pch,useCex=useCex,tit1=tit1,useCol=useCol,yLab=yLab,cexLa=cexLa,limM=limM,limA=limA,ptBg=ptBg)}  
    #check# plot(Avalue ~Mvalue, data=MAP14a$merg)
    #check# plot(MAP14a$merg[,"Avalue"], MAP14a$Mvalue)
    #check# plot(MAP14a$merg[,"Avalue"], MAP14a$Mvalue, pch=MAP14a$pch , cex=MAP14a$cex , main=MAP14a$tit1, ylab="A-value", col=MAP14a$useCol)
  } 


  if(length(Mvalue) >0) {
    pl1 <- try(graphics::plot(Mvalue ~ merg[,"Avalue"], pch=pch, cex=useCex, main=tit1,
      xlab="A-value", col=useCol, ylab=yLab, cex.lab=cexLa, xlim=limA,ylim=limM, pt.bg=ptBg), silent=TRUE)
    if(inherits(pl1, "try-error")) { Mvalue <- NULL; message("UNABLE TO PLOT !!  check plotting device ...")} }
  if(debug) {message(fxNa,"MAP14b"); MAP14b <-  list(merg=merg,Mvalue=Mvalue,filtFin=filtFin,passAll=passAll,pch=pch,useMar=useMar,FCthrs=FCthrs,useCol=useCol,useComp=useComp,useCompNa=useCompNa,subTxt=subTxt,singleCompSetup=singleCompSetup,inversedComp=inversedComp,setup=setup )} 

  if(length(Mvalue) >0) {                        # useCompNa
    sTxt <- if(length(subTxt) ==1) subTxt else { if(!singleCompSetup) { if(length(names(useComp)) >0) paste0(names(useComp),"; ",collapse="") else paste0(if(inversedComp) "inv.","comp=",useComp,"; ")}}
    sTxt <- paste0(sTxt,"n=",length(Mvalue),
      if(!all(is.na(c(FCthrs)))) paste(";",sum(passAll, na.rm=TRUE),"(color) points passing",
        if(!is.na(FCthrs)) paste0("(FCthr=", as.character(FCthrs),")") ))
    if(!singleCompSetup && "pwGrpNa" %in% names(setup) && length(dim(setup$pwGrpNa))==2) { sTxt <- paste(rownames(setup$pwGrpNa)[useComp],";", sTxt) }

    graphics::mtext(sTxt,cex=0.75,line=0.2)
    if(!all(is.na(c(FCthrs)))) {
      if(debug) message(fxNa," n=",length(Mvalue),"  FCthrs=",as.character(FCthrs),"  filt.ini=", sum(filtFin, na.rm=TRUE),
        "  passAll=",sum(passAll,na.rm=TRUE)," ; range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=TRUE),3))," ;  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
      ## previous versions : adjusted abline by adding  diff(graphics::par("usr")[1:2])/1000)  
      graphics::abline(v=c(-1,1)* log2(FCthrs) , col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
    if(debug) {message(fxNa,"MAP15"); MAP15 <-  list(merg=merg,Mvalue=Mvalue,filtFin=filtFin,FCthrs=FCthrs,annotColumn=annotColumn,yLab=yLab,tit1=tit1,col=col,useCol=useCol,assAll=passAll,grayIncrem=grayIncrem,cexPt=cexPt,annColor=annColor,limM=limM,limA=limA,ptBg=ptBg,cexLa=cexLa,passAll=passAll,namesNBest=namesNBest)}
#    if(sum(passFdr, na.rm=TRUE) >0) {
#      if(plotAsFdr) {
#        if(any(passAll, na.rm=TRUE)) graphics::abline(h=-1*log10(max(merg[passAll,"FDR"], na.rm=TRUE)) -diff(graphics::par("usr")[3:4])/400, col=grDevices::rgb(0.87,0.72,0.72), lty=2)
#      } else {
#        if(plotAsFdr) graphics::mtext("Note, that FDR and p-value may not correlate perfectly, thus points may appear at good p-value but finally don't get retained",line=5,cex=0.65, side=1)
#        pRa <- range(merg[which(passFdr),"Avalue"], na.rm=TRUE)
#        graphics::abline(h=pRa[1] +diff(graphics::par("usr")[3:4])/400, col=grDevices::rgb(0.87,0.72,0.72), lty=2) }}

    ## add names to best points
    if(length(namesNBest) >0) namesNBest <- wrMisc::naOmit(namesNBest)
    if(length(namesNBest) >0) {
      if(any(sapply( c("pass","passThr","signif","passFC"), identical, namesNBest))) namesNBest <- sum(passAll)
      if(!is.integer(namesNBest)) namesNBest <- try(as.integer(namesNBest), silent=TRUE)
      if(inherits(namesNBest, "try-error") || isTRUE(any(is.na(namesNBest)))) { namesNBest <- NULL
        message(fxNa,"Unable to understand argument 'namesNBest', must be integer or 'pass','passThr' or 'signif' (for display of labels to points in figure)") }
    }
    if(debug) {message(fxNa,"MAP16"); MAP16 <-  list(merg=merg,Mvalue=Mvalue,filtFin=filtFin,FCthrs=FCthrs,annotColumn=annotColumn,yLab=yLab,tit1=tit1,col=col,useCol=useCol,assAll=passAll,grayIncrem=grayIncrem,cexPt=cexPt,annColor=annColor,limM=limM,limA=limA,ptBg=ptBg,cexLa=cexLa,passAll=passAll,namesNBest=namesNBest)}
    if(length(namesNBest) >0) {
      if(namesNBest >0 && any(passAll, na.rm=T)) {
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
        if(length(noNa) >0 && all(annotColumn %in% colnames(merg))) names(tmP)[useL2[noNa]] <- merg[useLi[useL2[noNa]], wrMisc::naOmit(match(annotColumn[-1], colnames(merg)))[1]]
        if(length(NbestCol) <1) NbestCol <- 1
        if(is.null(names(tmP[useL2]))) {if(!silent) message(fxNa,"No names available for displaying names of best in plot")
        } else { displTx <- wrMisc::naOmit(match(annotColumn[2:4], colnames(merg)))
          displTx <- if(length(displTx) >0) cbind(merg[useLi[useL2], displTx], rownames(merg)[useLi[useL2]]) else as.matrix(useLi[useL2])
          chNa <- rowSums(!is.na(displTx)) >0
          if(any(chNa)) displTx[which(chNa),1] <- "unknown"
          displTx <- apply(displTx,1, function(x) wrMisc::naOmit(x)[1])
          graphics::text(Mvalue[useLi[useL2]] +xOffs, yOffs + merg[useLi[useL2],"Avalue"],
            names(tmP)[useL2], cex=cexTxLab, col=NbestCol, adj=0) }
      }
    }
    if(debug) {message(fxNa,"MAP17"); MAP17 <-  list(pch=pch,passAll=passAll,merg=merg,annot=annot,annotColumn=annotColumn)}

    ## legend (if multiple symbols)
    pch[which(is.na(pch))] <- -2
    ch1 <- unique(pch)
    if(length(ch1) >1 && sum(passAll) >0) {
      legInd <- which(!duplicated(merg[which(passAll), annotColumn[1]], fromLast=FALSE))
      if(length(legInd) <1 && !silent) message(fxNa,"Trouble ahead : Can't find non-duplicated ",annotColumn[1]," for ",sum(passAll)," points passing thresholds ! (ie as 'legInd')")
      legPch <- pch[which(passAll)[legInd]]
      legCol <- useCol[which(passAll)[legInd]]
      legBg <- ptBg[which(passAll)[legInd]]
      if(alph2 <1) {legCol <- substr(legCol,1,7); legBg <- substr(legBg,1,7)}  # reset to no transparency
      legLab <- merg[which(passAll)[legInd], annotColumn[1]]
      chNa <- is.na(legLab)
      if(any(chNa)) legLab[chNa] <- "NA"
      legOr <- if(length(legLab) >1) order(legLab) else 1   # not used so far
      legLoc <- checkForLegLoc(cbind(Mvalue, Avalue), sampleGrp=legLab, showLegend=FALSE)
      legCex <- stats::median(c(useCex, cexTxLab, 1.2), na.rm=TRUE)
      if(length(legLab) >0) { chLeg <- try(graphics::legend(legLoc$loc, legend=legLab, col=legCol, text.col=1, pch=legPch,
        if(length(ptBg) >0) pt.bg=ptBg, cex=legCex, pt.cex=1.2*legCex, xjust=0.5, yjust=0.5), silent=TRUE)  # as points
        if(inherits(chLeg, "try-error") && !silent) message(fxNa,"Note: Failed to add legend .. ",chLeg) }
    } else legCol <- NULL
    if(debug) {message(fxNa,"MAP18")}
    ## arrow for expected ratio
      ## how to extract in smartest way ??

    drawArrow <- if(length(expFCarrow) >0) {isTRUE(as.logical(expFCarrow[1])) || grepl("^[[:digit:]]+$", expFCarrow[1])} else FALSE

    if(drawArrow) {
      regStr <-"[[:space:]]*[[:alpha:]]+[[:punct:]]*[[:alpha:]]*"
      if(isTRUE(expFCarrow)) {
        ## automatic extraction of FC
        if(debug) message(fxNa," .. ",names(useComp),"  -> ",paste(unlist(strsplit(names(useComp), "-")), collapse=" "),"\n")
        expM <- sub(paste0("^",regStr),"", sub(paste0(regStr,"$"), "", unlist(strsplit(names(useComp), "-"))))    # assume '-' separator from pairwise comparison
        arrCol <- 1
        chN2 <- try(as.numeric(expM), silent=TRUE)
        if(inherits(chN2, "try-error")) { drawArrow <- FALSE #expM <- arrCol <- NA
        } else {
          expM <- log2(chN2[2] / chN2[1]); arrCol <- 1          # transform automatic extracted to ratio   # expFCarrow <- c(expM, 1)
        }
        msg <- "Argument 'expFCarrow' was set to TRUE; "
        if(!silent) message(fxNa, msg, if(drawArrow) "Extract concentration values of group-names for calculating ratios" else "Failed automatic extraction of FC-values (possibly no numeric values in setup)")
      } else {
        expM <- try(as.numeric(expFCarrow[1]), silent=TRUE)
        if(inherits(expM, "try-error")) {
          if(!silent) message(fxNa,"Argument 'expFCarrow' was given, but no numeric content in 1st place found,  can't draw FC-arrow")     #understand 1st value of'expFCarrow' (should be numeric-like)
          drawArrow <- FALSE
        } else arrCol <- if(length(expFCarrow) >1) expFCarrow[2] else 1
      }
      if(length(expFCarrow) >1) {arrCol <- expFCarrow[2]            # 2nd value of expFCarrow for color
      } else arrCol <- 1       # autom : use 3rd if possible
    }

    if(drawArrow) {
      if(is.finite(expM)) {
        figCo <- graphics::par("usr")                         #  c(x1, x2, y1, y2)
        figRa <- diff(range(figCo[1:2]))*0.1
        if(expM > figCo[2] +figRa || expM < figCo[1] -figRa) {
          if(!silent) message(fxNa,"Can't draw arrow, ",round(expM,2)," is too far outside the plotting frame")
          expM <- NA }
      }
      if(is.finite(expM)) {
        arr <- c(0.019,0.14)                                  # start- and end-points of arrow (as relative to entire plot)
        graphics::arrows(expM, figCo[3] + arr[1]*(figCo[4] -figCo[3]), expM, figCo[3] + arr[2]*(figCo[4] -figCo[3]),
          col=arrCol, lwd=1, length=0.1)
        graphics::mtext(paste("expect at",signif(expM,3)), at=expM, side=1, adj=0.5, col=arrCol, cex=cexLa*0.7, line=-0.9)
      } else { if(!silent) message(fxNa,"Unable to draw arrow for expexted M-value)") }
    }
    if(debug) {message(fxNa,"MAP19")}
    tmp <- try(graphics::par(mar=opar$mar, cex.main=opar$cex.main, las=opar$las), silent=TRUE)  # resetto previous

    ## export data used for plotting
    if(returnData) {
      merg1 <- try(cbind(merg, passThrs=passAll), silent=TRUE)
      if(!inherits(merg1, "try-error")) merg <- merg1
      merg <- merg[,-1*c(1, ncol(merg)-1)]        # remove col 'ID' 'redundant' & 'pch'
      annCo <- wrMisc::naOmit(match(annotColumn, colnames(merg)))
      if(length(annCo) >0) cbind(merg[,annCo],  merg[,-annCo]) else merg }
  } }
}
    
     
   