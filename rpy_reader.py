#==== PREAMBLE: python packages import ====
import sys
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages 
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.lib import grid

#==== PREAMBLE: r packages import ====
base=rpackages.importr("base")
ggplot2=rpackages.importr("ggplot2")
Cairo=rpackages.importr("Cairo")
igraph=rpackages.importr("igraph")
#Necessary to import graphical devices
grdevices=rpackages.importr('grDevices')
#R print function
rprint = robjects.globalenv.get("print")

grid.activate()

# !!! Import options??? !!!
#==== PREAMBLE ENDS ====

#---- Importing functions ----
definitions=base.source("definitions.R")

#---- Is the first argument defined? ----

if len(sys.argv)<2:
    print("No arguments given")
    PrintHelp=robjects.r('PrintHelp')
    PrintHelp(str(0))
    sys.exit(-1)

if sys.argv[1] and sys.argv[1]=="99":
    ballgown=rpackages.importr("ballgown")
    dplyr=rpackages.importr("dplyr",on_conflict="warn")
    genefilter=rpackages.importr("genefilter")

    robjects.r(
        '''
        print("Importing data")                                                                                                                           
        ##args[2] Mapping file                                                                                                                            
        ##args[3] Folder with data                                                                                                                        
        ##args[4] Common Part                                                                                                                             
        
        ##phenodata: reports name of folders where input are stored and the tissue they refer to.                                                         
        phenodata<-read.csv("phenodata.csv","header"=TRUE)                                                                                                        
                                                                                                                                                      
        #ballgown class: imports data from the data from outside folders                                                                                  
        data <- ballgown( dataDir = "/mnt/disk2/bufalo/data/SAMPLES",
        samplePattern = "ERR",                                                                                                             
        bamfiles = NULL,                                                                                                                     
        pData = phenodata,                                                                                                                   
        verbose = TRUE,                                                                                                                      
        meas = "all")                                                                                                                        
        #Issue: pdata in ballgow costructor doesn't work,                                                                                                 
        #Proceding in making it "manually"                                                                                                                
        pData(data)=phenodata
        #Getting only the genes with a significat variance between expressions in tissues                                                                 
        data.fil = subset(data,"rowVars(texpr(data))>1",genomesubset=TRUE)                                                                                
    
        save(data.fil,file="data.fil.RData")                                                                                                              
    
        ''')
    print("Your data are stored in data.fil.RData")


#---- Loading Data ----
try:
    print("Loading transcripts data")
    base.load("data.fil.RData")
    print("done")
except:
    print("ERROR: No data.fil.RData available: run this script with argument 99")
    exit -3

#---- Sanity Check ----
robjects.r(
    '''
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
    '''
)

#==== FUNCTIONS RUN ====

transcripts=robjects.r('transcripts')
phenodata=robjects.r('phenodata')

if str(sys.argv[1])=="1":
    if len(sys.argv)==3:
        Plotter=robjects.r('Plotter')
        Plot=Plotter(sys.argv[2],transcripts)
        rprint(Plot)
        
if sys.argv[1]=="2":
    if len(sys.argv)==4:
        SearchByTissue=robjects.r('SearchByTissue')
        Tissue=SearchByTissue(str(sys.argv[2]),str(sys.argv[3]),transcripts)
    else:
        SearchByTissue=robjects.r('SearchByTissue')
        Tissue=SearchByTissue(str(sys.argv[2]),"",transcripts)

    out_file=open("SearchTissue.tsv","w",5000000)
    out_file.write(str(Tissue))
    out_file.close
        
if sys.argv[1]=="3":
    SearchByGene=robjects.r('SearchByGene')
    Gene=SearchByGene(str(sys.argv[2]),transcripts)

    out_file=open("SearchGene.tsv","w",5000000)
    out_file.write(str(Gene))
    out_file.close

if sys.argv[1]=="4":
    SearchByFeature=robjects.r('SearchByFeature')
    data_fil=robjects.r('data.fil')
    if len(sys.argv)==4:
        Feature=SearchByFeature(str(sys.argv[2]),str(sys.argv[3]),data_fil,phenodata)
    else:
        Feature=SearchByFeature(str(sys.argv[2]),'',data_fil,phenodata)


    out_file=open("SearchFeature.tsv","w",5000000)
    out_file.write(str(Feature))
    out_file.close

if sys.argv[1]=="5":
    SearchByDiffFoldExpr=robjects.r('SearchByDiffFoldExpr')
    if len(sys.argv)==5:
        Fold=SearchByDiffFoldExpr(str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),transcripts)
    else:
        Fold=SearchByDiffFoldExpr(str(sys.argv[2]),str(sys.argv[3]),"",transcripts)
    
    out_file=open("SearchFold.tsv","w",5000000)
    out_file.write(str(Fold))
    out_file.close

if sys.argv[1]=="6":
    SearchTranscriptGroup=robjects.r('SearchTranscriptGroup')
    Module=SearchTranscriptGroup(str(sys.argv[2]),str(sys.argv[3]),transcripts)

    out_file=open("SearchModule.tsv","w",5000000)
    out_file.write(str(Module))
    out_file.close

if sys.argv[1]=="7":
    Network=robjects.r('Network')
    NET=Network(str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),transcripts,str(sys.argv[5]),str(sys.argv[6]))
    rprint(NET)
