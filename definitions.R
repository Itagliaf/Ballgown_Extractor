library(ballgown)	
library(dplyr)
library(genefilter)

lotter2<-function(GENE, MEAS="FPKM", SAMPLES=pData(bg)$ids,BASEDIR=getwd(), bg)
{
    ##plots all transcripts relative to a single gene relatively to all tissues.
    
    ##vvvv add random string from stringi package?vvvv
    file_hash<-paste(format(Sys.time(), "%d%H%M%s"))
    ##^^^^ add random string from stringi package?^^^^
    
    
    gene<-toupper(GENE)
    
    ## vvvv MARCONI CAN'T CREATE PNG, switching to pdf
    
    out_file1 <- paste(GENE, MEAS, file_hash,"png", sep=".")
    out_file2 <- paste(BASEDIR,out_file1,sep="/")
    png(out_file2,width=1080,height=720)
    
    par(mar = rep(2, 4))    
    plotTranscripts(GENE,bg, meas=MEAS, sample=SAMPLES, labelTranscripts=T, legend=T, colorby="transcript")
    
    dev.off()
    
    return(0)
}




SearchByCondition<-function(COMBO_CONDITIONS, GENE=NULL, bg)
    
    ##this function results in a table containing the FPKM and cov values for a gene in a specific subset of samples. If no gene name is provided, it prints the values for all the genes in a specific subset of samples or in a specific sample. 
    
{
    bg_subset <-subset(bg, COMBO_CONDITIONS, genomesubset=FALSE)
    
    Transcripts<-bg_subset@expr$trans
    GENE<-toupper(GENE)
    
    if (!is.null(GENE))
    {
        if(GENE %in% Transcripts$gene_id | GENE %in% Transcripts$gene_name)
        {
            print(paste("Displaying FKPM and cov values for gene", GENE, "in the sample:", pData(bg_subset)$ids))
            attach(Transcripts)
            final <- Transcripts[which(gene_id==GENE | gene_name==GENE),]
            detach(Transcripts)
        }
        else
        {   
            print(paste("Displaying FKPM and cov values for Samples", pData(bg_subset)$ids))
            return(Transcripts)
        }  
    }
    return(final)
}

SearchByGene<-function(GENE,bg)
    ##subsets the dataframe to extract only lines with a certain gene_name (gene symbol) or gene_id
    
{
    print(paste("Search By Gene", GENE))
    Transcripts<-bg@expr$trans
    if(GENE %in% Transcripts$gene_id | GENE %in% Transcripts$gene_name)
    {
        id <- unique(Transcripts[which(Transcripts$gene_id==GENE | Transcripts$gene_name==GENE),9])
        
        ##vvvv NOT actually needed for genes but its necessary for lowercase
        ##Name<-toupper(Name)
        ##Gene<-tolower(GENE)
        final <- gexpr(bg)[id, ]
        return(final)
    } 
    else 
    { 
        print("Invalid Gene Name inserted")
    }
    
}

SearchByTranscript<-function(TRANSCRIPT_ID, bg)
    ##subsets the dataframe to extract only lines with a certain gene_name (gene symbol)
{
    print(paste("Search By Transcript", TRANSCRIPT_ID))
    total <-bg@expr$trans[,-c(1:5,7:8)]
    rownames(total)<-total$t_name
    final<-total[TRANSCRIPT_ID,]
    return(final)
}


SearchGeneIsoforms<-function(GENE,bg)
    ##Given a gene id this function prints the expression values for its isoforms (when present). Furthermore, it adds metainformation regarding the genomic coordinates and the nucleotide length.
{
    print(paste("Search Isoforms for gene", GENE))
    total <-bg@expr$trans
    final<-total[(total$gene_id%in%GENE | total$gene_name %in% GENE),-c(1, 7, 8)]
    
    return(final)
}



SearchByFeature<-function(GENE,FEATURE,bg)                
    ##subsets the dataframe to extract only lines with a certain gene_name and a certain gene feature (exon or intron)
{
    print("Search By Gene Feature")
    
    print(GENE)
    print(FEATURE)
    
    phenodata<-pData(bg)
    
    ##the code from search by gene name is repeated 2 times and merged
    
    db<-bg@expr
    results=NULL                    
    if(GENE %in% db$trans$gene_name)            
    {
        total <-db$trans
        results <- total[total$gene_name%in%GENE,]
        ## I extract the chromosome where is located the gene and the coordinates of the considered gene
        chr_position<-results$chr
        gene_begin<-min(results$start)
        gene_end<-max(results$end)
    }
    else             
    {                
        print("Not Found")
        result=NaN
    }
    
    print("====STARTS SECOND IF====")
    
    if(FEATURE=="exon")
    {
        final<-db$exon[which(db$exon$chr==chr_position & db$exon$start >= gene_begin & db$exon$end <= gene_end),]
    }
    else
    {
        final<-db$intron[which(db$intron$chr==chr_position & db$intron$start >= gene_begin & db$intron$end <= gene_end),]
    }
    
    print("==== FINISH SECOND IF====")
    ## Change col names
    exclude<-grep("mrcount",colnames(final))
    
    final<-final[,-exclude]            
    
    exclude<-grep("ucount",colnames(final))
    
    final<-final[,-exclude]            
    
    NamesOriginal<-colnames(final)
    
    Sample_Names <- phenodata[,2]
    
    print("==== STARTS FOR LOOP ====")
    for(i in c(0:nrow(phenodata)))            
    {
        
        Folders=phenodata[i,1]                
        rcount="rcount"
        Sample_RC=paste(rcount,Folders,"sep"=".")
        if(Sample_RC%in%colnames(final))
        {
            pos<-match(Sample_RC,NamesOriginal)
            sample<- phenodata[i,2]
            str_sample <- paste("Read Count",toString(sample),sep=" ")
            
            NamesOriginal[pos]<- str_sample
        }
    }
    colnames(final)<-NamesOriginal
    return(final)
}


SearchByDiffFoldExpr<-function(COMBO_CONDITIONS, COVARIATE, FEATURE, bg)
{
    ##ipotesi di due o più variabili che arrivano già in questo formato dall'utente per il subsetting (ad esempio COMBO_CONDITIONS='time_h==336 & treatment=="Irradiated"'):
    
    bg_subset <-subset(bg,COMBO_CONDITIONS, genomesubset=FALSE)
    bg_filt<-subset(bg_subset,"rowVars(texpr(bg))>1",genomesubset=TRUE)
    
    ##Saving some metainformation:
    details_trans <-bg_filt@expr$trans[,c(1:10)]
    details_exon <-bg_filt@expr$exon[,c(1:5)] 
    details_intron <- bg_filt@expr$intron[,c(1:5)]
    
    
    if(FEATURE=="transcript") 
    {
        stats <- stattest(bg_filt, feature=FEATURE, meas='FPKM', covariate=COVARIATE, getFC=T)
        ##adding metainfo and sorting by qvalue:
        final_stats <- arrange(merge(stats, details_trans, by.x=2, by.y=1, all.x=T, all.y=F), qval)
        
        ##rearranging columns to show:
        return(final_stats[,c(10,14,3:9,11,12)])
    }
    else if(FEATURE=="gene") 
    {    
        stats <- stattest(bg_filt, feature=FEATURE, meas='FPKM', covariate=COVARIATE, getFC=T)
        
        final <- arrange(stats,qval)
        return(final[,c(2:5)])
    }
    else if(FEATURE=="exon")
    {
        stats <- stattest(bg_filt, feature=FEATURE, covariate=COVARIATE, getFC=T)
        final <- arrange(merge(stats,details_exon, by.x=2, by.y=1, all.x=T, all.y=F),qval)
        return(final[,c(2:9)])
    }
    else if(FEATURE=="intron")
    {
        stats <- stattest(bg_filt, feature=FEATURE, meas='rcount', covariate=COVARIATE, getFC=T)
        final <- arrange(merge(stats,details_intron, by.x=2, by.y=1, all.x=T, all.y=F),qval)
        return(final[,c(2:9)])
    }
}

StatsFiltering<-function(STATS, Q_THRESHOLD=0.05, P_THRESHOLD=0.05, MIN_FOLD_CHANGE=2)
{
    ##Remove NAs:                                                                                                                                         
    stats_without_na <- STATS[!is.na(STATS$pval), ]
    
    ##Print some summary stats:                                                                                                                           
    p.val_filtering <- sum(stats_without_na$pval< P_THRESHOLD)
    q.val_filtering <- sum(stats_without_na$qval< Q_THRESHOLD)
    fc_filtering <- sum(stats_without_na$fc> log2(MIN_FOLD_CHANGE + 1) | stats_without_na$fc< log2((1/MIN_FOLD_CHANGE)+1))
    combination_of_filters <- sum(stats_without_na$pval< P_THRESHOLD & stats_without_na$qval< Q_THRESHOLD & stats_without_na$fc> log2(MIN_FOLD_CHANGE + 1) | stats_without_na$fc< log2((1/MIN_FOLD_CHANGE)+1))
    
    pval_string<-paste("There are",p.val_filtering,"entries satisfying the desired p-value threshold",sep=" " )
    qval_string<-paste("There are",q.val_filtering,"entries satisfying the desired q-value threshold",sep=" " )
    fc_string<-paste("There are",fc_filtering,"entries satisfying the desired fold change threshold",sep=" " )
    combination_filters_string<-paste("There are", combination_of_filters, "entries satisfying the desired combination of fold change, p-value and q-value thresholds",sep=" " )
    
    print(pval_string)
    print(qval_string)
    print(fc_string)
    print(combination_filters_string)
    
    ##filter and show the data frame with the specified thresholds of interest:                                                                           
    filt_stats_table <- stats_without_na[(stats_without_na$pval< P_THRESHOLD & stats_without_na$qval< Q_THRESHOLD & (stats_without_na$fc > log2(MIN_FOLD_CHANGE +1) | stats_without_na$fc <log2(1/MIN_FOLD_CHANGE)+1)) ,]
    
    return(filt_stats_table)
}

getGenes<-function(bg)
{
    names<-geneNames(bg)
    ids <- geneIDs(bg)
    genes <- unique(cbind(names,ids))
    return(genes)
}

getTranscrip<-function(bg)
{
    names<-transcriptNames(bg)
    return(names)
}


Gene_Plotter_By_Group<-function(GENE, MEAS="FPKM", GROUPVAR=NULL, BASEDIR=getwd(),bg)
{
    ##plots all transcripts relative to a single gene grouped by a covariate (a colname of pData(bg)).
    
    if(is.null(GROUPVAR))
    {
        stop("Covariate variable is not defined")
    }
    
    ##vvvv add random string from stringi package?vvvv
    file_hash<-paste(format(Sys.time(), "%d%H%M%s"))
    ##^^^^ add random string from stringi package?^^^^
    
    gene<-toupper(GENE)
    out_file_1 <- paste(GENE, MEAS, GROUPVAR, file_hash,"png", sep=".")
    out_file_2 <- paste(BASEDIR,out_file_1,sep="/")
    
    png(out_file_1,width=1080,height=720)
    
    par(mar = rep(2, 4))    
    plotMeans(gene,bg, meas=MEAS, groupvar=GROUPVAR)
    
    dev.off()    
    return(0)
}


SearchTranscriptGroup<-function(NAME,MODULESFILE,bg)
{
    ##function returns the genes belonging to the same expression module as the gene of interest.     
    
    ##MODULESFILE: file containig informations about coexpression                                                                                     
    ##modules (TOM files from WGCNA pipeline)                                                                                                         
    ##It's the file OK_Gene_Modules.Rdata obtainde with coexpression.R                                                                                
    
    load(MODULESFILE)
    
    transcripts<-bg@expr$trans
    
    ##vvvvfind the position of the gene in the Expression Data                                                                                        
    
    pos<-match(NAME,colnames(datExpr))
    geneModule<-moduleLabels[pos]
    
    ## datExpr is an object included in MODULESFILE and contains the                                                                                  
    ## expression information                                                                                                                         
    ##^^^^find the position of the gene in the Expression Data                                                                                        
    
    ##extracting all genes belonging to the module                                                                                                    
    Gene_Module<-names(datExpr)[moduleLabels==geneModule]
    subsetted<-transcripts[transcripts$gene_id %in% Gene_Module,]
    
    return(subsetted[,c(1:10)])
}

Network<-function(QUERY,MODULES,corr=0.9,results=50,bg)
{
    load(bg)
    Transcripts=bg@expr$trans

    ##file hash for output
    File_Hash<-paste(format(Sys.time(), "%d%H%M"),"png",sep=".")

    if(file.exists(MODULES))
    {
        load(MODULES)
    }else{
        print("File not found. Makes sure that the file containing the modules is correct")
    }

    ##== END ==

    QUERY_position<-match(QUERY,colnames(datExpr))
    QUERY_module<-moduleLabels[QUERY_position]


    module_line<-which(moduleLabels %in% QUERY_module)
    datExpr_module<-datExpr[,module_line]

    Corr_Genes<-cor(datExpr_module)
    print("Correlation matrix created")
    colnames(Coor_Genes)

    diag(Coor_Genes)<-0
    Coor_Genes[abs(Coor_Genes)<corr]=0

    row_col=match(QUERY,colnames(Coor_Genes))
    cor_vec=abs(Coor_Genes[,row_col])
    C<-names(sort(abs(Coor_Genes[,row_col]),decreasing=TRUE)[c(0:results)])

    if(QUERY%in%C){
        print("QUERY gene is in C")
    }else{
        C[results]<-QUERY
    }
    D<-Coor_Genes[C,C]

    ##D/40 to stress distances between nodes
    D<-D/40

    print("Create graph")

    if(QUERY_type=="gene_name")
    {

        graph<-graph_from_adjacency_matrix(D,mode="max",weighted=TRUE,diag=FALSE)

        V(graph)$gene=as.character(transcripts[match(C,transcripts$gene_id),]$gene_name)
        #V(graph)$gene=as.character(transcripts[match(C,transcripts$t_name),]$gene_name)
        coul=rainbow(n=length(V(graph)$gene))
        my_color=coul[as.numeric(as.factor(V(graph)$gene))]

        V(graph)$color<-"white"
        V(graph)[QUERY]$color<-"red"
        V(graph)$shape<-"vrectangle"
        graph<-simplify(graph)

        ##plotting

        out_file=paste("Network",QUERY,File_Hash,sep="_")

        ##png(out_file,width=1080,heigh=720)

        plot(graph, vertex.shape="vrectangle", vertex.label.color=my_color,label.degree="-p/2")

        cols=2
        print("Fixing Legend")
        if(length(C)>70)
        {
            cols=3

        }else{
            cols=2
        }
        legend(x=-2,y=0.8,legend=levels(as.factor(V(graph)$gene)),
               col = coul ,
               bty = "n",
               pch=20 ,
               pt.cex = 3,
               cex = 1,
               text.col=coul,
               horiz = FALSE,
               ncol=cols,
               inset = c(0.1, 0.1))
        ##dev.off()
    }

    return(graph)
}


