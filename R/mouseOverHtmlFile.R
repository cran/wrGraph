#' Create mouse-over interactive html-pages (with links) 
#'
#' @description
#' This function allows generating html pages with interactive mouse-over to display information for the points of the plot and www-links when clicking based on embedded png file.  
#' Basically, an html page will be generated which contains a call to display to an image file specified in \code{pngFileNa} and in the body below pixel-coordinated will be 
#' given for disply of information at mouse-over and embedded links.
#' 
#' @details
#' Basically theer are two options for defining the path to the image embedded : 
#' 1) Absolute path : I turn you can moove the html to different locations, as long as it still can see the png-file the image can be displayed. However, this may 
#'  not be any more the case when the html file is sent to another person. If the png-file is accessible as url, it should be easily visible. 
#' 2) Relative path : The simplest case would be to give only the file-name with no path at all, thus the png-file is supposed to be in the same directory as the html-file. 
#' This option is very 'transportable'. 
#' Basically the same applies to the clickable links which may be provided. In high-throughput biology one typically points here to data-bases accessible
#' over the internet where urls to specific pages. With UniProt such links can easily be constructed when using protein identifiers as rownames.
#'	 
#' @param myCoor (matrix or data.frame) with initial x&y coordinates of points for plot; with IDs (1st column !!) & coordinates (2nd & 3rd col), data for mouse-over & link (4th & 5th); 
#'  NOTE : if 'colNa' NOT given, colnames of 'myCoor' will be inspected & filtered (columns of non-conform names may get lost) !!!  
#'  Associated with (already existing) figure file 'pngFileNa' and make html page where points may be indicated by mouse-over  
#' @param pngFileNa (character, length=1) filename for complementary png figure (must already exist)
#' @param HtmFileNa (character, length=1) filename for html file produced 
#' @param mouseOverTxt (character, length=1) text for interactive mouse-over in html, if \code{NULL}, will use col specified by 1st 'colNa' or (if NULL) rownames of 'myCoor'
#' @param displSi (integer, length=2)  size of image ('pngFileNa') at display in html (width,height), see also \code{\link[graphics]{par}} 
#' @param colNa (character) if not \code{NULL} min length of 3 to custom specify the column-names to be used : 1st for mouse-over and 2nd+3rd for coordinates associated (and optional 4th for links) 
#' @param tit (character) title to be displayed on top of figure 
#' @param myHtmTit (character) title of Html page; 'htmlExt' .. checking and correcting filename-extension (only main Html page) 
#' @param myComment (character) modify comment embedded in html-document 
#' @param textAtStart (character) text in html before figure 
#' @param textAtEnd (character) text in html after figure 
#' @param pxDiam (integer, length=1) diameter for mouse-over tip to appear (single val or vector), simpler version/solution than with 'Tooltip' package 
#' @param addLinks (character) for clickable links, either 1) vector of links or 2) single character-chain to be used for pasting to rownames (eg http://www.uniprot.org/uniprot/)  
#'  or 3) \code{TRUE} to check presence of 4th name specified in 'colNa' to be useed as columname from 'myCoor'   dominates over eventual presence of 4th name in 'colNa' 
#' @param linkExt (character) if specified : links will get specified ending, define as \code{NULL} or "" for taking 'addLinks' asIs  
#' @param htmlExt (character, length=1) extension used when making html files  
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @return plot
#' @seealso \code{\link{convertPlotCoordPix}}; use \code{\link[wrMisc]{htmlSpecCharConv}} to convert special characters for proper html display 
#' @examples
#' ## Note, this example writes files to R's tempdir,
#' ## Otherwise, if you simply work in the current directory without spcifying paths you'll 
#' ##  get an html with relatove paths, which simply needs the png file in the same path 
#' df1 <- data.frame(id=letters[1:10],x=1:10,y=rep(5,10),mou=paste("point",letters[1:10]),
#'   link=file.path(tempdir(),paste(LETTERS[1:10],".html",sep="")),stringsAsFactors=FALSE)  
#' ## here we'll use R's tempdir, later you may want to choose other locations
#' pngFile <- file.path(tempdir(),"test01.png")
#' png(pngFile,width=800, height=600,res=72)
#' ## here we'll just plot a set of horiontal points ...
#' plot(df1[,2:3],las=1,main="test01")
#' dev.off()
#' ## Note : Special characters should be converted for display in html pages during mouse-over
#' library(wrMisc)
#' df1$mou <- htmlSpecCharConv(df1$mou)
#' ## Let's add the x- and y-coordiates of the points in pixels to the data.frame
#' df1 <- cbind(df1,convertPlotCoordPix(x=df1[,2],y=df1[,3],plotD=c(800,600),plotRes=72))
#' head(df1)
#' ## Now make the html-page allowing to display mouse-over to the png made before
#' htmFile <- file.path(tempdir(),"test01.html")
#' mouseOverHtmlFile(df1,pngFile,HtmFileNa=htmFile,pxDiam=15,
#'   textAtStart="Points in the figure are interactive to mouse-over ...",
#'   textAtEnd="and/or may contain links")
#' ## We still need to make some toy links
#' for(i in 1:nrow(df1)) cat(paste("point no ",i," : ",df1[i,1]," x=",df1[i,2]," y=",
#'   df1[i,3],sep=""), file=df1$link[i]) 
#' ## Now we are ready to open the html file using any browser
#' \dontrun{ 
#' browseURL(htmFile)
#' }
#' @export
mouseOverHtmlFile <- function(myCoor, pngFileNa, HtmFileNa=NULL, mouseOverTxt=NULL, displSi=c(800,600),
	colNa=NULL, tit="", myHtmTit="", myComment=NULL, textAtStart=NULL, textAtEnd=NULL, pxDiam=5,
    addLinks=NULL, linkExt=NULL, htmlExt="htm", callFrom=NULL, silent=FALSE){
  ## make html file/output where (www-)links & supplemental information is accessible at mouse-over on image/plot (eg IDs/names in plot)
  ##  assume, that initial png is already made and coordinates are already converted to pixel level of png-image
  ##  'displSi'  size of image ('pngFileNa') at display in html (width,height)
  ## 'addLinks' for clickable links, either 1) vector of links or 2) single char-chain to be used for pasting to rownames (eg http://www.uniprot.org/uniprot/)
  ##   or 3) TRUE to check presence of 4th name specified in 'colNa' to be useed as columname from 'myCoor'
  ##   dominates over eventual presence of 4th name in 'colNa'
  myCoorTy <- NULL
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="mouseOverHtmlFile")
  if(length(dim(myCoor)) !=2) stop(" Expecting matrix or data.frame")
  if(nrow(myCoor) <1) stop(" 'myCoor' seems to be empty !")
  if(!is.data.frame(myCoor)) myCoor <- as.data.frame(myCoor,stringsAsFactors=FALSE)
  if(is.null(myComment)) myComment <- c(" Produced by R using createHtmlWithPointsIdentif()", " from package WRmisc, ",Sys.Date())
  .corPath <- function(x,asHtml=TRUE) {
    ## correct mixed slash and backslash in file path
    if(length(grep("ming.32",R.Version()$platform)) >0) {
      x <- gsub("\\\\","/",x)  #"
      if(asHtml & length(grep("[[:upper:]]:",substr(x,1,2))) >0) {
        x <- paste("file:///",x,sep="") }
    } else if(asHtml & length(grep("^/",x)) >0) x <- paste("file:///",x,sep="")
    x }      
  ## check colNa
  colN2 <- rep(NA,5)  
  names(colN2) <- c("ID","x","y","mouseOver","link")
  potIDname <- wrMisc::.plusLowerCaps(c("ID","Id","Ident","Identifier","UniqID","uniqID","UniqueID"))
  potXname <- wrMisc::.plusLowerCaps(c("xPix","xPred","coorX","htmlX","xCoor","dataX","X"))
  potYname <- wrMisc::.plusLowerCaps(c("yPix","yPred","coorY","htmlY","yCoor","dataY","Y"))
  potMouseOvname <- wrMisc::.plusLowerCaps(c("mouseOver","mouseInfo","Mouse","Mou","Hint","Name","fullName","FullName","Combined","CustName","custName","Info"))
  potLinkname <- wrMisc::.plusLowerCaps(c("LINK","Link","Li","Http","WWW","AddLink","addLink"))
  ## determine column nmaes to be extracted
  if(length(colNa) ==2) colN2[2:3] <- colNa else if(length(colNa) >2) colN2[1:length(colNa)] <- colNa
  if(length(colNa) <1) {
    if(ncol(myCoor)==2) colN2[2:3] <- colnames(myCoor) else {
      remColNa <- colnames(myCoor)
      if(any(potIDname %in% remColNa)) {aa <- .serachColName(remColNa,potIDname,plusLowerCaps=FALSE,returnList=TRUE)
        remColNa <- aa$remainNa; colN2[1] <- aa$foundNa}
      if(any(potXname %in% remColNa)) {aa <- .serachColName(remColNa,potXname,plusLowerCaps=FALSE,returnList=TRUE)
        remColNa <- aa$remainNa; colN2[2] <- aa$foundNa}
      if(any(potYname %in% remColNa)) {aa <- .serachColName(remColNa,potYname,plusLowerCaps=FALSE,returnList=TRUE)
        remColNa <- aa$remainNa; colN2[3] <- aa$foundNa}
      if(any(potMouseOvname %in% remColNa)) {aa <- .serachColName(remColNa,potMouseOvname,plusLowerCaps=FALSE,returnList=TRUE)
        remColNa <- aa$remainNa; colN2[4] <- aa$foundNa}
      if(any(potLinkname %in% remColNa)) {aa <- .serachColName(remColNa,potLinkname,plusLowerCaps=FALSE,returnList=TRUE)
        remColNa <- aa$remainNa; colN2[5] <- aa$foundNa}}}
  ## select data to be used
  myCoor <- data.frame(myCoor[,colN2[which(!is.na(colN2))]],stringsAsFactors=FALSE)              # set proper order of cols
  if(is.na(colN2[1])) {
    myCoor$ID <- if(is.null(rownames(myCoor))) 1:nrow(myCoor) else rownames(myCoor)
    colN2[1] <- "ID"
    if(!silent) message(fxNa," using column '",colnames(myCoor[colN2[1]]),"' for mouse-over")
  }
  if(!identical(mouseOverTxt,FALSE)) if(is.null(mouseOverTxt)) {
    if(is.na(colN2[4])) {
      myCoor$mouseOver <- myCoor[,colN2[1]]                                             # nothing specified, use names as default
      if(!silent) message(fxNa," using column '",colnames(myCoor[colN2[1]]),"' for mouse-over")
      colN2[4] <- "mouseOver" }
  } else { if(length(mouseOverTxt) ==nrow(myCoor) & length(unique(mouseOverTxt)) > 1) {
    myCoor$mouseOver <- mouseOverTxt } else {
      if(!silent) message(fxNa,"Ignoring invalid entry for 'mouseOverTxt' (expecting length ",nrow(myCoor)," but found ",length(mouseOverTxt),")")
      myCoor$mouseOver <- myCoor[,"ID"]
    } }
  ##
  if(identical(addLinks,TRUE)) {            # colNa has priority over IDs
    myCoor$link <- myCoor[,if(is.na(colN2[5])) colN2[1] else colN2[5]]
    colN2[5] <- "addLinks" } else {
   if(length(addLinks) >0) {if(length(addLinks) ==nrow(myCoor) & length(unique(addLinks)) > 1 &max(nchar(addLinks),na.rm=TRUE) >0) {
     myCoor$link <- addLinks
     colN2[5] <- "addLinks" } else {
     if(!silent) message(fxNa," invalid entry for 'addLinks' (expecting length ",nrow(myCoor)," but found ",length(addLinks),")")}}} # no default
  ## check extensions of specified links (ie myCoor$link) : if 'linkExt' specified (& longer than 0 char) add to this ending if not present
  if(length(linkExt) >0) if(nchar(linkExt) >0) {
    chExt <- grep(paste(linkExt,"$",sep=""),myCoor$link)
    if(length(chExt) < nrow(myCoor) & length(chExt) >0) myCoor$link[chExt] <- paste(myCoor$link[chExt],linkExt,sep="")
  }
  ## prepare for making html file
  if(!file.exists(pngFileNa)) stop("Cannot find file which should be used for embedding image into html !")
  msg <- " 'displSi' : Expecting numeric vector of lengt 2 (for display size in px in html)  !"
  if(!is.numeric(displSi) | length(displSi) <2) stop(msg)
  if(length(HtmFileNa) !=1) HtmFileNa <- pngFileNa
  baseFiNa <- sub(".PNG$","",sub(".png$","",sub(".htm$","",sub(".html$","",HtmFileNa))))
  if(nchar(HtmFileNa)== nchar(baseFiNa)) {
    htmlExt <- if(length(htmlExt) <0) "" else htmlExt[1]                     # allow file wo extesion if empty argument 'htmlExt'
    if(!silent) message(fxNa," setting file-name + extension to : ",baseFiNa,".",htmlExt)
  } else {
    htmlExt <- substr(HtmFileNa,unlist(regexec("\\.htm",HtmFileNa)),nchar(HtmFileNa))}
  HtmFileNa <- wrMisc::.checkFileNameExtensions(baseFiNa,htmlExt)       # check file extensions for HtmFileNa & pngFileNa
  .convTxtToHtmPar <- function(txt){    # convert character vector to paragraphs  <p>My paragraph.</p>
    txt <- as.character(txt)
    nLi <- length(txt)
    apply(matrix(c(rep("<p>",nLi), txt,rep("</p>",nLi)),nrow=nLi),1,paste,collapse="") }
  ## main, ie html creation
  htmVec <- c('<!DOCTYPE html>','<html lang="en">','<head>','<meta charset="utf-8">')
  htmTit <- paste(c('<title>',myHtmTit,'</title>'),collapse="")
  htmVec <- c(htmVec,htmTit,"</head>","<body>")
  htmCom <- paste(c("<!-- ",myComment,"-->"),collapse="")
  htmGraTit <- if(is.null(tit)) NULL else paste(c("<h2>",tit,"</h2>"),collapse="")           # graphic title
  htmVec <- c(htmVec,htmCom,htmGraTit)
  if(!is.null(textAtStart)) htmVec <- c(htmVec,.convTxtToHtmPar(textAtStart))
  htmImg <- paste(c('<img src="',.corPath(pngFileNa),'" alt="wrGraph_imageForMouseOver" usemap="#colormap" style="width:',
    displSi[1],'px;height:',displSi[2],'px">'),collapse="")
  htmVec <- c(htmVec,htmImg,'<map name="colormap">')
  ar1 <- '<area title="'
  ar3 <- 'shape="circle" coords="'
  ar5 <- ' alt="'
  ar7 <- '"'
  htmCor <- data.frame(ar1,na1=myCoor[,colN2[4]],'" ',ar3,corX=myCoor[,colN2[2]],
    ',',corY=myCoor[,colN2[3]],',',diam=pxDiam,naZ='"',stringsAsFactors=FALSE)
  if(!is.na(colN2[5])) htmCor <- data.frame(htmCor[,-1*ncol(htmCor)],na2='" href="',na3=.corPath(myCoor$link),'"',stringsAsFactors=FALSE)
  htmCor <- cbind(htmCor,last=" >")
  htmCor <- as.character(apply(htmCor,1,paste,collapse=""))
  htmVec <- c(htmVec,htmCor,"</map>")
  if(!is.null(textAtEnd)) htmVec <- c(htmVec,.convTxtToHtmPar(textAtEnd))
  htmVec <- c(htmVec,"</body>","</html>")
  if(is.null(HtmFileNa)) HtmFileNa <- paste(sub(".png$","",pngFileNa),".html",sep="")
  if(file.exists(HtmFileNa) & !silent) message(fxNa," BEWARE, file '",HtmFileNa,"' will be overwritten !")
  tryWrite <- try(cat(paste(htmVec,collpse="\n"),file=HtmFileNa))
  if(class(tryWrite) =="try-error") warning(fxNa," PROBLEM : couldn't write Html file '",
    HtmFileNa,"' ! (file open ?  check path,rights etc)")
}        

#' @export
.serachColName <- function(x,searchColNa, plusLowerCaps=TRUE,returnList=TRUE,callFrom=NULL) {
  ## 'x' character vector of column-names to inspect, or matrix/data.frame where colnames will be extracted/inspected
  ## 'searchColNa' (character) 
  fxNa <- wrMisc::.composeCallName(callFrom,newNa=".serachColName")
  x <- if(length(dim(x)) >1)  colnames(x) else as.character(x)
  if(length(x) <1) return(NULL) else {
  errMsg <- "argument 'searchColNa' is emty or all NA !"
  if(length(searchColNa) <1) stop(errMsg)
  chNa <- is.na(searchColNa)
  if(any(chNa)) {if(all(chNa)) stop(errMsg) else searchColNa <- searchColNa[which(!chNa)]}
  if(plusLowerCaps) searchColNa <-  wrMisc::.plusLowerCaps(searchColNa)
  out <- wrMisc::naOmit(match(searchColNa,x))
  if(length(out) <1) stop("none of the terms found")
  if(returnList) list(foundNa=x[out[1]],remainNa=x[-out]) else out[1] }}
  
