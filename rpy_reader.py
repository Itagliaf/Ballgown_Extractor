#==== PREAMBLE: python packages import ====
import sys
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages 

#==== PREAMBLE: r packages import ====
base=rpackages.importr("base")
ggplot2=rpackages.importr("ggplot2")
Cairo=rpackages.importr("Cairo")
igraph=rpackages.importr("igraph")

# !!! Import options??? !!!
#==== PREMABLE ENDS ====

#---- Importing functions ----
base.source("definitions.R")

#---- creating data.fil.Rdata
#print(sys.argv[1])

# if sys.argv[1]=="99":
#     ballgown=rpackages.importr("ballgown")
#     dplyr=rpackages.importr("dplyr",on_conflict="warn")
#     genefilter=rpackages.importr("genefilter")

#     robjects.r(
#         '''
#         print("Importing data")                                                                                                                           
#         ##args[2] Mapping file                                                                                                                            
#         ##args[3] Folder with data                                                                                                                        
#         ##args[4] Common Part                                                                                                                             
        
#         ##phenodata: reports name of folders where input are stored and the tissue they refer to.                                                         
#         phenodata<-read.csv("phenodata.csv","header"=TRUE)                                                                                                        
                                                                                                                                                      
#         #ballgown class: imports data from the data from outside folders                                                                                  
#         data <- ballgown( dataDir = "/home/itagliaferri/Documents/Progeti/Chillemi/SAMPLES",
#         samplePattern = "ERR",                                                                                                             
#         bamfiles = NULL,                                                                                                                     
#         pData = phenodata,                                                                                                                   
#         verbose = TRUE,                                                                                                                      
#         meas = "all")                                                                                                                        
#         #Issue: pdata in ballgow costructor doesn't work,                                                                                                 
#         #Proceding in making it "manually"                                                                                                                
#         pData(data)=phenodata
#         #Getting only the genes with a significat variance between expressions in tissues                                                                 
#         data.fil = subset(data,"rowVars(texpr(data))>1",genomesubset=TRUE)                                                                                
    
#         save(data.fil,file="data.fil.RData")                                                                                                              
    
#         ''')
#     print("Your data are stored in data.fil.RData")


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

options = {0 : PrintHelp(0),
           1 : Plotter(sys.argv[2],)
           }

options[sys.argv[1]]()
