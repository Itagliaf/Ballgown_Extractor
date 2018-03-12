#==== PREAMBLE: python packages import ====
import sys
import datetime
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
rplot = robjects.globalenv.get("plot")
rlegend = robjects.globalenv.get("legend")
grid.activate()

# !!! Import options??? !!!
#==== PREAMBLE ENDS ====

#==== USEFULL FUNCTIONS ====

def exp_spacer_up():
    """
    To highlight a wrong called function (head)
    """
    print("\n=====vvvvvvvvvv=====")

def exp_spacer_down():
    """
    To highlight a wrong called function (bottom)
    """
    print("=====^^^^^^^^^^=====\n")
#==== END =====


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

#to format output files name
now = datetime.datetime.now()
time=now.strftime("%Y%m%d%H%M")


# if str(sys.argv[1])=="1":
#     print("\n Function: Plotter")
#     if len(sys.argv)==3:
#         Plotter=robjects.r('Plotter')
#         Plot=Plotter(sys.argv[2],transcripts)
        
#     else:
#         print("\n ATTENTION: Parsing Arguments Error")
#         exp_spacer_up()
#         print("Argument 1: gene symbol of gene you want to analyze or file containing the gene names to be analyzed (sep='\\n')")
#         exp_spacer_down()
#         sys.exit("Exiting")

#     out_name="Plotter_"+sys.argv[2]+"_"+time+".png"
#     grdevices.png(file=out_name, width=1080, height=720)
#     rprint(Plot)
#     grdevices.dev_off()
    
if str(sys.argv[1])=="1":
    print("\n Function: Plotter")
    if len(sys.argv)==3:
        Plotter=robjects.r('Plotter2')
        Plot=Plotter2(sys.argv[2],transcripts)    
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: gene symbol of gene you want to analyze or file containing the gene names to be analyzed (sep='\\n')")
        exp_spacer_down()
        sys.exit("Exiting")
        
if sys.argv[1]=="2":
    print("\n Function: SearchByTissue")
    if len(sys.argv)==4:
        SearchByTissue=robjects.r('SearchByTissue')
        Tissue=SearchByTissue(str(sys.argv[2]),str(sys.argv[3]),transcripts)
    elif len(sys.argv)==3: 
        SearchByTissue=robjects.r('SearchByTissue')
        Tissue=SearchByTissue(str(sys.argv[2]),"",transcripts)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 2: tissue name to be analyzed")
        print("Argument 3: (optional): gene symbol of gene to be analyzed")
        exp_spacer_down()
        sys.exit("Exiting")

    out_name="SearchTissue_"+sys.argv[2]+"_"+time+".tsv"
    out_file=open(out_name,"w",5000000)
    out_file.write(str(Tissue))
    out_file.close
        
if sys.argv[1]=="3":
    print("\n Function: SearchByGene")
    if len(sys.argv)==3:
        SearchByGene=robjects.r('SearchByGene')
        Gene=SearchByGene(str(sys.argv[2]),transcripts)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 2: gene symbol to be analyzed")
        exp_spacer_down()
        sys.exit("Exiting")

    out_name="SearchGene_"+sys.argv[2]+"_"+time+".tsv"
    out_file=open(out_name,"w",5000000)
    out_file.write(str(Gene))
    out_file.close

if sys.argv[1]=="4":
    print("\n Function: SearchByFeature")
    SearchByFeature=robjects.r('SearchByFeature')
    data_fil=robjects.r('data.fil')
    if len(sys.argv)==4:
        Feature=SearchByFeature(str(sys.argv[2]),str(sys.argv[3]),data_fil,phenodata)
    elif len(sys.argv)==3:
        Feature=SearchByFeature(str(sys.argv[2]),'',data_fil,phenodata)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: Gene symbol to be analyzed")
        print("Argument 2 (Optional): Gene feature to be analyzed (exon or intron)")
        exp_spacer_down()
        sys.exit("Exiting")
        
    out_name="SearchFeature_"+sys.argv[2]+"_"+sys.argv[3]+"_"+time+".tsv"
    out_file=open(out_name,"w",5000000)
    out_file.write(str(Feature))
    out_file.close

if sys.argv[1]=="5":
    print("\n Function: SearchByDiffFoldExpr")
    SearchByDiffFoldExpr=robjects.r('SearchByDiffFoldExpr')
    if len(sys.argv)==5:
        Fold=SearchByDiffFoldExpr(str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),transcripts)
    elif len(sys.argv)==4:
        Fold=SearchByDiffFoldExpr(str(sys.argv[2]),str(sys.argv[3]),"",transcripts)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: First Tissue to be analyzed")
        print("Argument 2: Second Tissue to be analyzed")
        print("Argument 3 (optional): a particular gene to be analyzed (gene symbol)")
        exp_spacer_down()
        sys.exit("Exiting")
        
    out_name="SearchFold_"+sys.argv[2]+"_"+sys.argv[3]+"_"+time+".tsv"
    out_file=open(out_name,"w",5000000)
    out_file.write(str(Fold))
    out_file.close

if sys.argv[1]=="6":
    print("\n Function: SearchTranscritpGroup")
    if len(sys.argv)==5:
        SearchTranscriptGroup=robjects.r('SearchTranscriptGroup')
        Module=SearchTranscriptGroup(str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),transcripts)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: ")
        print("Argument 2: Gene name (eg 'gene77') or gene symbol (eg 'STK33') to be analyzed")
        print("Argument 3: File containing informations about coexpression (gene modules)")
        print("Argument 4: TOM files from WGCNA pipeline")
        exp_spacer_down()
        sys.exit("Exiting")
        
    out_name="SearchExpressionModule_"+sys.argv[2]+"_"+time+".tsv"
    out_file=open(out_name,"w",5000000)
    out_file.write(str(Module))
    out_file.close

if sys.argv[1]=="7":
    print("\n Function: Network")
    if len(sys.argv)==7:
        Network=robjects.r('Network')
        NET=Network(str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),transcripts,str(sys.argv[5]),str(sys.argv[6]))
        rplot(NET)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: gene name or gene symbol is used? Possible values: symbol or gene_name")
        print("Argument 2: Gene name (eg 'gene77') or gene symbol (eg 'STK33') to be analyzed")
        print("Argument 3: File containing informations about coexpression (gene modules)")
        print("""Argument 4: correlation threashold value to show a connection.
NB: this value is considered ad absolute value: both positive and negative connections are considered""")
        print("Argument 5: number of genes (nodes) to be mantained")
        exp_spacer_down()
        sys.exit("Exiting")


if sys.argv[1]=="8":
    print("\n Function: Gene Fold Change Between Tissues")
    if len(sys.argv)==3:
        GeneFoldTissue=robjects.r("GeneFoldTissue")
        GeneFold=GeneFoldTissue(argv[2],transcritps)
    else:
        print("\n ATTENTION: Parsing Arguments Error")
        exp_spacer_up()
        print("Argument 1: gene symbol (Eg stk33, abcd1) to be searched")
        exp_spacer_down()
        sys.exit("Exiting")
