library(WGCNA)
library(ballgown)

##!!!!!!!!!!!!!!!!!!!!!!

##the parts between

###vvvvvvvvvvvvvvvvvvv
##and
###^^^^^^^^^^^^^^^^^^^

##should be changed according to your system and necessity

## this "hardcoding problem" resolved in the next releases
##!!!!!!!!!!!!!!!!!!!!!!


## The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

## Threads number: a thread per core is the best choiche

###vvvvvvvvvvvvvvvvvvv
tn<-20
###^^^^^^^^^^^^^^^^^^^
enableWGCNAThreads(nThreads =20)


## For names
source("definitions.R")

load("data.fil.RData")


transcripts <-data.fil@expr$trans
phenodata<-read.csv("phenodata.csv","header"=TRUE)
colnames(transcripts)<-NameFormatter(transcripts,phenodata)

##vvv Now is GENE expression
## HARCODED must change
ok<-gexpr(data.fil)

##transpose to put gene expression as variable and tissue as US
##!! datExpr0 is the name used in the Tutorial
datExpr0=as.data.frame(t(ok))

gsg = goodSamplesGenes(datExpr0, verbose = 3);

##if this is TRUE ok, else remove "bad genes"
gsg$allOK

if (!gsg$allOK)
{
    ## Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    ## Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

collectGarbage()

##plot to see if there's outliar

sampleTree = hclust(dist(datExpr0), method = "average");
## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
## The user should change the dimensions if the window is too large or too small.

##opens a window, dimensions in inches (WHHHHYYYYY???!?!?!?!?)
## sizeGrWindow(12,9)

## par(cex = 0.6);
## par(mar = c(0,4,2,0))
## plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
##      cex.axis = 1.5, cex.main = 2)

## Plot a line to show the cut (here 100000)
## abline(h = 100000, col = "red");

## Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 10000, minSize = 10)
table(clust)

## clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

## Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

## Plot the results:
## sizeGrWindow(9, 5)
## par(mfrow = c(1,2));
## cex1 = 0.9;
##                                         # Scale-free topology fit index as a function of the soft-thresholding power
## plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
##      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
##      main = paste("Scale independence"));
## text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
##      labels=powers,cex=cex1,col="red");
##                                         # this line corresponds to using an R^2 cut-off of h
## abline(h=0.90,col="red")
##                                         # Mean connectivity as a function of the soft-thresholding power
## plot(sft$fitIndices[,1], sft$fitIndices[,5],
##      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
##      main = paste("Mean connectivity"))
## text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


## Construction the gene network and identifing modules
net = blockwiseModules(datExpr, power = 7,maxBlockSize=30000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM_Complete",
                       verbose = 3)

moduleLabels = net$colors
save(moduleLabels,datExpr,file = "Gene_Modules.Rdata")
