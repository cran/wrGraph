## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----setup2, echo=TRUE, eval=FALSE--------------------------------------------
# if(!requireNamespace("wrMisc", quietly=TRUE)) install.packages("wrMisc")   # required underneath
# install.packages("wrGraph")

## ----setup3, echo=TRUE, eval=FALSE--------------------------------------------
# if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")   # also required
# if(!requireNamespace("FactoMineR", quietly=TRUE)) install.packages("FactoMineR")
# if(!requireNamespace("factoextra", quietly=TRUE)) install.packages("factoextra")

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE------------------------
suppressPackageStartupMessages({
    library(wrMisc)
    library(wrGraph)
    library(FactoMineR)
    library(factoextra)
})

## ----vignBrowse, echo=TRUE, eval=FALSE----------------------------------------
# # access vignette :
# browseVignettes("wrGraph")    #  ... and the select the html output

## ----setup4, echo=TRUE--------------------------------------------------------
library("wrMisc")
library("wrGraph")
# This is version no:
packageVersion("wrGraph")

## ----partitionPlot1, out.width="110%", out.heigth="110%", echo=TRUE-----------
## as the last column of the Iris-data is not numeric we choose -1
(part <- partitionPlot(ncol(iris)-1))
layout(part)
for(i in colnames(iris)[-5]) hist(iris[,i], main=i)

## ----LegLoc1, echo=TRUE-------------------------------------------------------
dat1 <- matrix(c(1:5,1,1:5,5), ncol=2)
grp <- c("abc","efghijk")
(legLoc <- checkForLegLoc(dat1, grp)) 

# now with more graphical parameters (using just the best location information)
plot(dat1, cex=2.5, col=rep(2:3,3),pch=rep(2:3,3))
legLoc <- checkForLegLoc(dat1, grp, showLegend=FALSE)
legend(legLoc$loc, legend=grp, text.col=2:3, pch=rep(2:3), cex=0.8)

## ----Hist1, out.width="110%", out.heigth="110%", echo=TRUE--------------------
set.seed(2016); dat1 <- round(c(rnorm(200,6,0.5), rlnorm(300,2,0.5), rnorm(100,17)),2)
dat1 <- dat1[which(dat1 <50 & dat1 > 0.2)]
histW(dat1, br="FD", isLog=FALSE)

## ----Hist2, out.width="110%", out.heigth="110%", echo=TRUE--------------------
## view as log, but x-scale in linear
histW(log2(dat1), br="FD", isLog=TRUE, silent=TRUE)

## ----Hist4, out.width="110%", out.heigth="150%", echo=TRUE--------------------
## quick overview of distributions  
layout(partitionPlot(4))
for(i in 1:4) histW(iris[,i], isLog=FALSE, tit=colnames(iris)[i])

## ----Hist5, out.width="110%", out.heigth="150%", echo=TRUE--------------------
layout(1)
plot(iris[,1:2], col=rgb(0.4,0.4,0.4,0.3), pch=16, main="Iris Data")
legendHist(iris[,1], loc="br", legTit=colnames(iris)[1], cex=0.5)
legendHist(iris[,2], loc="tl", legTit=colnames(iris)[2], cex=0.5)

## ----vioplot1, fig.height=6, fig.width=8, echo=TRUE---------------------------
vioplotW(iris[,-5], tit="Iris-Data Violin Plot")

## ----vioplot2, fig.height=6, fig.width=8, echo=TRUE---------------------------
## less smoothing
vioplotW(iris[,-5], tit="Iris-Data Violin Plot ('nervous')", hh=0.15)

## ----vioplot3, fig.height=6, fig.width=6, echo=TRUE---------------------------
## less smoothing
vioplotW(iris[,-5], tit="Paired Iris-Data Violin Plot ", halfViolin="pairwise")

## ----cumFrqPlot1, echo=TRUE---------------------------------------------------
cumFrqPlot(iris[,1:4])

## ----imageW, echo=TRUE--------------------------------------------------------
par(mar=c(4, 5.5, 4, 1))  
imageW(as.matrix(iris[1:40,1:4]), transp=FALSE, tit="Iris-Data (head)")

## ----imageW2, fig.height=2.5, fig.width=9, fig.align="center", echo=TRUE------
imageW(as.matrix(iris[1:20,1:4]), latticeVersion=TRUE, transp=FALSE, col=c("blue","red"), 
  rotXlab=45, yLab="Observation no", tit="Iris-Data (head)")

## ----imageW3, fig.height=6, fig.width=5, echo=TRUE----------------------------
ma1 <- matrix(-7:16,nc=4,dimnames=list(letters[1:6],LETTERS[1:4]))
ma1[1,2:3] <- 0
ma1[3,3] <- ma1[3:4,4] <- NA
   
imageW(ma1, latticeVersion=TRUE, col=c("blue","grey","red"), NAcol="grey92", 
  rotXlab=0, cexDispl=0.8, tit="Balanced color gradient")   

## ----imageW4, fig.height=6, fig.width=5, echo=TRUE----------------------------
imageW(ma1, latticeVersion=TRUE, col=c("blue","grey","red"), NAcol="grey92", 
  rotXlab=0, nColor=8, cexDispl=0.8, tit="Balanced color gradient")   

## ----cumulCountPlot1, echo=TRUE-----------------------------------------------
thr <- seq(min(iris[,1:4]), max(iris[,1:4])+0.1,length.out=100)
irisC <- sapply(thr,function(x) colSums(iris[,1:4] < x))
irisC <- cbind(thr,t(irisC))

  head(irisC)
staggerdCountsPlot(irisC[,], countsCol=colnames(iris)[1:4], tit="Iris-Data")
staggerdCountsPlot(irisC[,], varCountNa="Sepal", tit="Iris-Data")
staggerdCountsPlot(irisC[,], varCountNa="Sepal", tit="Iris-Data (log-scale)", logScale=TRUE)

## ----plotBy2Groups1, echo=TRUE------------------------------------------------
dat <- iris[which(iris$Species %in% c("setosa","versicolor")),]
plotBy2Groups(dat$Sepal.Length, gl(2,50,labels=c("setosa","versicolor")),
  gl(20,5), yLab="Sepal.Length")

## ----plotLinReg1, echo=TRUE---------------------------------------------------
plotLinReg(iris$Sepal.Length, iris$Petal.Width, tit="Iris-Data")

## ----PCA1, fig.height=7, fig.width=7, echo=TRUE-------------------------------
## the basic way
iris.prc <- prcomp(iris[,1:4], scale.=TRUE)
biplot(iris.prc)              # traditional plot

## ----PCA3a, echo=TRUE---------------------------------------------------------
## via FactoMineR
chPa <- c(requireNamespace("FactoMineR", quietly=TRUE), requireNamespace("dplyr", quietly=TRUE), 
  requireNamespace("factoextra", quietly=TRUE), requireNamespace("ggpubr", quietly=TRUE) )

## ----PCA3b, echo=TRUE---------------------------------------------------------
if(all(chPa)) {
  library(FactoMineR); library(dplyr); library(factoextra)
  iris.Fac <- PCA(iris[,1:4],scale.unit=TRUE, graph=FALSE)
  iris.Fac2 <- try(fviz_pca_ind(iris.Fac, geom.ind="point", col.ind=iris$Species, palette=c(2,4,3), 
    addEllipses=TRUE, legend.title="Groups"))
  if(inherits(iris.Fac2, "try-error")) message("Problem running factoextra ..") else {
    iris.Fac2 <- try(plot(iris.Fac2))  
    if(inherits(iris.Fac2, "try-error")) message("A problem occured when trying to print from factoextra !!") }    
} else message("You need to install packages 'dplyr', 'FactoMineR' and 'factoextra' for this figure ! ",
  "Check the messages of library() if there are other packages that may be missing (and need to get installed first)")  


## ----PCA4, echo=TRUE----------------------------------------------------------
## via wrGraph, similar to FactoMineR but with bagplots
plotPCAw(t(as.matrix(iris[,-5])), gl(3,50,labels=c("setosa","versicolor","virginica")),
  tit="Iris Data", rowTyName="types of leaves", suplFig=FALSE, cexTxt=1.3, rotatePC=2)

## ----PCA5, fig.height=12, fig.width=9, fig.align="center", echo=TRUE----------
## including 3rd component and Screeplot
plotPCAw(t(as.matrix(iris[,-5])), gl(3,50,labels=c("setosa","versicolor","virginica")),
  tit="Iris Data PCA", rowTyName="types of leaves", cexTxt=2)


## ----PCA6, fig.height=7.5, fig.width=9, fig.align="center", echo=TRUE---------
## creat copy of data and add rownames
irisD <- as.matrix(iris[,-5])
rownames(irisD) <- paste(iris$Species, rep(1:50,3), sep="_")
## plot
plotPCAw(t(irisD), gl(3,50,labels=c("setosa","versicolor","virginica")), tit="Iris Data PCA", 
  rowTyName="types of leaves", suplFig=FALSE, cexTxt=1.6, rotatePC=2, pointLabelPar=list(textCex=0.48))

## ----MA0, fig.height=6, fig.width=8, fig.align="center", echo=TRUE------------
## toy data
set.seed(2005); mat <- matrix(round(runif(2400),3), ncol=6)
mat[11:90,4:6] <- mat[11:90,4:6] +round(abs(rnorm(80)),3)
mat[11:90,] <- mat[11:90,] +0.3
dimnames(mat) <- list(paste("li",1:nrow(mat),sep="_"),paste(rep(letters[1:2],each=3),1:6,sep=""))
## assume 2 groups with 3 samples each
matMeans <- round(cbind(A=rowMeans(mat[,1:3]), B=rowMeans(mat[,4:6])),4)

## ----MA1, fig.height=6, fig.width=8, fig.align="center", echo=TRUE------------
## now we are ready to plot, M-values can be obtained by subtracting thr group-means
MAplotW(M=matMeans[,2] -matMeans[,1], A=rowMeans(mat)) 

## ----MA4, echo=TRUE-----------------------------------------------------------
## assume 2 groups with 3 samples each and run moderated t-test (from package 'limma')
tRes <- wrMisc::moderTest2grp(mat, gl(2,3), addResults=c("FDR","Mval","means"))

## ----MA5,fig.height=6, fig.width=8, fig.align="center",  echo=TRUE------------
## convenient way, change fold-change threshold to 2x and mark who is beyond :
MAplotW(tRes, FCth=2, namesNBest="passFC")    

## ----Volc1, fig.height=6, fig.width=8, fig.align="center", echo=TRUE----------
## let's generate some toy data
set.seed(2005); mat <- matrix(round(runif(900),2), ncol=9)
rownames(mat) <- paste0(rep(letters[1:25],each=4), rep(letters[2:26],4))
mat[1:50,4:6] <- mat[1:50,4:6] + rep(c(-1,1)*0.1,25)
mat[3:7,4:9] <- mat[3:7,4:9] + 0.7
mat[11:15,1:6] <- mat[11:15,1:6] - 0.7

## assume 2 groups with 3 samples each
gr3 <- gl(3,3,labels=c("C","A","B"))
tRes2 <- moderTest2grp(mat[,1:6], gl(2,3))

VolcanoPlotW(tRes2)

# now with thresholds, labels and arrow for expected ratio
VolcanoPlotW(tRes2, FCth=1.3, FdrThrs=0.2, namesNBest="pass", expFCarrow=c(0.75,2))

## ----Volc2,  fig.height=6, fig.width=9.5, fig.align="center", echo=TRUE-------
## assume 3 groups with 3 samples each
tRes <- moderTestXgrp(mat, gr3)

layout(matrix(1:2, nrow=1))
VolcanoPlotW(tRes, FCth=1.3, FdrThrs=0.2, useComp=2)
VolcanoPlotW(tRes, FCth=1.3, FdrThrs=0.2, useComp=3)

## ----createHtmlWithPointsIdentif1, echo=TRUE----------------------------------
## Let's make some toy data 
df1 <- data.frame(id=letters[1:10], x=1:10, y=rep(5,10) ,mou=paste("point",letters[1:10]),
  link=file.path(tempdir(),paste(LETTERS[1:10],".html",sep="")),stringsAsFactors=FALSE)  
## here we'll use R's tempdir, later you may want to choose other locations
pngFile <- file.path(tempdir(),"test01.png")
png(pngFile, width=800, height=600, res=72)
## here we'll just plot a set of horiontal points ...
plot(df1[,2:3], las=1, main="test01")
dev.off()
## Note : Special characters should be converted for proper display in html during mouse-over
df1$mou <- htmlSpecCharConv(df1$mou)
## Let's add the x- and y-coordiates of the points in pixels to the data.frame
df1 <- cbind(df1, convertPlotCoordPix(x=df1[,2], y=df1[,3], plotD=c(800,600), plotRes=72))
head(df1)
## Now make the html-page allowing to display mouse-over to the png made before
htmFile <- file.path(tempdir(),"test01.html")
mouseOverHtmlFile(df1, pngFile, HtmFileNa=htmFile, pxDiam=15,
  textAtStart="Points in the figure are interactive to mouse-over ...",
  textAtEnd="and/or may contain links")
## We still need to make some toy links
for(i in 1:nrow(df1)) cat(paste("point no ",i," : ",df1[i,1]," x=",df1[i,2]," y=",
  df1[i,3],sep=""), file=df1$link[i]) 
## Now we are ready to open the html file using any browser ..
#from within R# browseURL(htmFile) 

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

