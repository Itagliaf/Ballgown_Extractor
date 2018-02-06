##---- Preamble: Importing Libraries----
library(ggplot2)
library(Cairo)
#antialiasing
options(shiny.usecairo=T)
##---------------------------------------

##----- Importing function definitions ----
source("definitions.R")
#------------------------------------------

##---- Begin execution ----
#Feeding arguments to the script
args <- commandArgs(TRUE)


##-- Importing files --
#Sanity check: there are arguments?
if (length(args)==0)
{
    PrintHelp(0)
}

if (args[1]==99 && length(args)<4)
{
    print("Importing data")
    print("* WARNING: insufficient arguments")
    print("Argmiument 2: Mapping file (phenodata file)")
    print("Argument 3: Folder containing the data")
    print("Argument 4: Common part of sample folders' name")
    stop("Insufficient Arguments")
}else if (args[1]==99 && length(args)==4)
{

    print("Importing data")

    library(ballgown)
    library(dplyr)
    library(genefilter)

    ##args[2] Mapping file
    ##args[3] Folder with data
    ##args[4] Common Part

    ##phenodata: reports name of folders where input are stored and the tissue they refer to.
    phenodata<-read.csv(args[2],"header"=TRUE)

    #ballgown class: imports data from the data from outside folders
    data <- ballgown( dataDir = args[3],
                 samplePattern = args[4],
                 bamfiles = NULL,
                 pData = phenodata,
                 verbose = TRUE,
                 meas = "all")
    #Issue: pdata in ballgow costructor doesn't work,
    #Proceding in making it "manually"
    pData(data)=phenodata

    #---- Data preparing ----
    #Getting only the genes with a significat variance between expressions in tissues
    data.fil = subset(data,"rowVars(texpr(data))>1",genomesubset=TRUE)

    save(data.fil,file="data.fil.RData")

    print("Your data are stored in data.fil.RData")
    stop("All done: exiting")
}

## -- Loading data --
## if data is already imported
tryCatch(
{
    load("data.fil.RData")
},
warning=function(w)
{
    print("ERROR: No data.fil.RData available: run the script with the firs argument=99")
},
warning=function(e)
{

    stop("Cannot find data.fil")
}
)

##Extracting only the transcripts: we have fpkm normalized of each gene inside samples
transcripts <-data.fil@expr$trans

if (length(transcripts)==0)
{
    stop("ERROR: the transcripts dataframe is empty. Is data.fil.RData empty?")
}

##phenodata (for names)
phenodata<-read.csv("phenodata.csv","header"=TRUE)

##NameFormatter(transcripts,phenodata)
colnames(transcripts)<-NameFormatter(transcripts,phenodata)


# -- Actual Funcions run -- 

##Function run: for each function there's a sanity check

File_Hash<-paste(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),"tsv",sep=".")

switch(args[1],
       "1"=if(length(args)<2)
           {
              print("Plotter: plots FPKM vaules of one or more genes in all tissue considered")
              print("Argument 2: File containig the genes names to be analyzed (sep=\n)")
              print("Exiting")
              stop("Insufficient Arguments")
           }
           else
           {
               print("Plotter: plots FPKM vaules of one or more genes in all tissue considered")
               Genes<-scan(args[2],what="character")
               Plotter(Genes,transcripts)
           },
       "2"=if(length(args)<2)
           {
               print("Search by Tissue")
               print("Argument 2: tissue name to be analyzed")
               print("Argument 3 (otpional): a particular gene to analyzed")
               print("Exiting")
               stop("Insufficient Arguments")
           }
           else
           {
               print("Search by Tissue")
               Tissue_Done<-SearchByTissue(args[2],args[3],transcripts)
               out_file=paste("SearchTissue",args[2],File_Hash,sep="_")
               write.table(Tissue_Done, out_file,row.names=FALSE,col.names=TRUE)
           },
       "3"=if(length(args)<2)
           {
               print("Search by Gene")
               print("Argument 2: gene name to be analyzed")
               print("Exiting")
               stop("Insufficient Arguments")
           }
           else
           {
               print("Search by Gene")
               Gene_Done<-SearchByGene(args[2],transcripts)
               out_file=paste("SearchGene",args[2],File_Hash,sep="_")
               write.table(Gene_Done, out_file,row.names=FALSE,col.names=TRUE)
           },
       "4"=if(length(args)<3)
           {
               print("Search by gene feature")
               print("Argument 2: gene name to be analyzed")
               print("Argument 3: gene feature to be analyzed")
               print("Exiting")
               stop("Insufficient Arguments")
           }else
           {
               print("Search by gene feature")
               ##args[2]=gene name
               ##args[3]=feature name
               Feature_Done<-SearchByFeature(args[2],args[3],data.fil,phenodata)
               out_file=paste("SearchFeature",args[2],args[3],File_Hash,sep="_")
               write.table(Feature_Done, out_file,row.names=FALSE,col.names=TRUE)
           },
       "5"=if(length(args)<3)
           {
               print("Search by Differenctila fold")
               print("Argument 2: first tissue to be analyzed")
               print("Argument 3: second tissue to be analyzed")
               print("Argument 4 (optional): a particular gene to be analyzed (name)")
               print("Exiting")
               stop("Insufficient Arguments")
           }else
           {
               print("Search By Differential Fold")
               ##args[2]=tissue1
               ##args[3]=tissue2
               ##args[4]=gene name (optional) to subset
               Fold_Done <- SearchByDiffFoldExpr(args[2],args[3],args[4],transcripts)
               out_file=paste("SearchFold",args[2],args[3],File_Hash,sep="_")
               write.table(Fold_Done, out_file,row.names=FALSE,col.names=TRUE)
           },
       "6"={
           Gene_Module<-SearchTranscriptGroup(args[2],args[3],transcripts)
           print(Gene_Module)
           }
      ,
       PrintHelp()
       
)
