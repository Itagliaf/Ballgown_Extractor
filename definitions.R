NameFormatter<-function(Transcripts,phenodata)
    ##formats in a better way the names of the columns

{
    ##takes the names from the transcripts dataframe
    NamesOriginal<-colnames(Transcripts)

    ##takes the names of the tissue from phenodata
    Tissue_Names<-phenodata[,2]
    
    ##for each row in the phenotipic data
    for (i in c(0:nrow(phenodata)))
    {
        ##take the name of the ith sample and add ".cov" or ".FPMK"
        Folders=phenodata[i,1]  
        Fpkm="FPKM"
        Cov="cov"
        TissueCol=paste(Fpkm,Folders,"sep"=".")
        TissueCov=paste(Cov,Folders,"sep"=".")

        ##If TissueCol is found in the names of the transcripts, sobstitute the name
        if(TissueCol%in%colnames(transcripts))
        {
            pos<-match(TissueCol,NamesOriginal)
            tissue<-phenodata[i,2]
            str_tissue<-paste(Fpkm,toString(tissue),sep=".")
            NamesOriginal[pos]<-str_tissue
        }

        ##If TissueCov is found in the names of the transcripts, sobstitute the name	
        if(TissueCov%in%colnames(transcripts))
        {
            pos<-match(TissueCov,NamesOriginal)
            tissue<-phenodata[i,2]
            str_tissue<-paste(Cov,toString(tissue),sep=".")
            NamesOriginal[pos]<-str_tissue
        }
        

    }
    ##returns a vector with the names of the tissues
    return(NamesOriginal)
}

PrintHelp<-function()
{
    print ("Choose a function asfollow:")
    print ("1 = Plotter Function")
    print ("2 = Search by Tissue")
    print ("3 = Search by gene")
    print ("4 = Search by gene feature")
    print ("5 = Search by Differential Fold")
    print ("99 = import data from a ballgonw object")
    stop("No arguments given")
}

SearchByTissue<-function(Tissue,Name)
    ##subsets the dataframe to extract only columns reguarding a certain tissue
{
    
    ##create FPKM.tissue and OUT.tissue strings
    paste("Plotting FKPM for tissue ",Tissue)
    Fpkm="FPKM"
    Out="OUT"
    TissueCol=paste(Fpkm,Tissue,"sep"=".")
    ##valutate the value of TissueCol and search it in transcripts
    temp<-quote(TissueCol)

    ##subset transcript with only the fpkm column
    TransFpkm<-transcripts[ ,eval(temp)]
    
    ##Add coverage column in the same way as before
    Cov="cov"
    TissueCol=paste(Cov,Tissue,"sep"=".")
    temp<-quote(TissueCol)
    TransCov<-transcripts[ ,eval(temp)]

    ##create a dataframe that has the id's of the genes and the fpkm data
    final<-merge(transcripts$gene_name,TransFpkm,by=0,all=TRUE)

    ##add the coverage values in a column
    final$Coverage=TransCov

    ##rename the columns
    colnames(final)<-c("ID","Gene_name","Fpkm","Coverage")
    
    if (!is.null(Name) & !is.na(Name))
    {
        if(Name %in% final$Gene_name)
        {
            attach(final)
            final <- final[which(Gene_name==Name),]
            detach(final)
        }
        else
        {
            stop(sprintf("Can't find gene: %s",Name))
        }  
    }
        
    
    return(final)
}

SearchByGene<-function(Name,Transcripts)
    ##subsets the dataframe to extract only lines with a certain gene_name (gene symbol)
{       
    ##if the gene is oresente in the column gene_name
    ##extract the line with the gene and return it
    if(Name %in% Transcripts$gene_name)
    {
        attach(Transcripts)
        final <- Transcripts[which(gene_name==Name),]
        detach(Transcripts)
    }
    else
    {
        stop(sprintf("Can't find gene: %s",Name))
    }
    return(final)
}

SearchByFeature<-function(Name,Feature,data.fil,Phenodata)
    ##subsets the dataframe to extract only lines with a certain gene_name and a certain gene feature (exon or intron)
    ##data.fil must be input (presents informations about the intron, exon and so on.
{
    ##the code from search by gene name is repetead 2 times and merged
    db<-data.fil@expr
    results=NULL
    if(Name %in% db$trans$gene_name)
    {

        transcripts <-db$trans
        attach(transcripts)
        results <- transcripts[which(gene_name==Name),]
        detach(transcripts)

        ## I extract the chromosome where is located the gene
        ##and the start-end position of the gene considered 
        chr_position<-results$chr
        gene_begin<-min(results$start)
        gene_end<-max(results$end)
    }
    else
    {
        print("Not Found")
        result=NaN
    }

    if(Feature=="exon")
    {
        final<-db$exon[which(db$exon$chr==chr_position & db$exon$start >= gene_begin & db$exon$end <= gene_end),]
    }
    else
    {
        final<-db$intron[which(db$intron$chr==chr_position & db$intron$start >= gene_begin & db$intron$end <= gene_end),]
    }

    ## Change name of the columns
    ## Other columns to be maintained?
    
    exclude<-grep("mrcount",colnames(final))
    final<-final[,-exclude]
    exclude<-grep("ucount",colnames(final))
    final<-final[,-exclude]

    print(head(final))

    NamesOriginal<-colnames(final)
    
    Tissue_Names <- phenodata[,2]
    for (i in c(0:nrow(phenodata)))
    {
        Folders=phenodata[i,1]
        rcount="rcount"
        Tissue_RC=paste(rcount,Folders,"sep"=".")
        if(Tissue_RC%in%colnames(final))
        {
            pos<-match(Tissue_RC,NamesOriginal)
            tissue <- phenodata[i,2]
            str_tissue <- paste("Read Count",toString(tissue),sep=" ")
            NamesOriginal[pos] <- str_tissue
        }
    }
    colnames(final)<-NamesOriginal
    return(final)
}

SearchByDiffFoldExpr<-function(Tissue1,Tissue2,Name,Transcripts)
    ##subsets the dataframe to extract only 2 tissues and confront their expression (fold expression)
{       
    ##Extract FPKM values from the tissues 
    Fpkm="FPKM"
    Out="OUT"
    TissueCol1=paste(Fpkm,Tissue1,"sep"=".")
    TissueCol2=paste(Fpkm,Tissue2,"sep"=".")

    temp1<-quote(TissueCol1)
    temp2<-quote(TissueCol2)
    
    TransFpkm1<-Transcripts[,eval(temp1)]
    TransFpkm2<-Transcripts[,eval(temp2)]

    ##Add covariance column

    Cov="cov"
    TissueCol1=paste(Cov,Tissue1,sep=".")
    TissueCol2=paste(Cov,Tissue2,sep=".")

    temp1<-quote(TissueCol1)
    temp2<-quote(TissueCol2)
    
    TransCov1<-Transcripts[,eval(temp1)]
    TransCov2<-Transcripts[,eval(temp2)]

    ##Creating a data table containing only the columns needed
    LittleData<-merge(Transcripts$gene_name,TransFpkm1,by=0,all=TRUE)
    LittleData$Transcript2=TransFpkm2

    colnames(LittleData)<-c("ID","Gene_name",Tissue1,Tissue2)

    if (!is.null(Name) & !is.na(Name))
    {
        if(Name %in% LittleData$Gene_name)
        {
            attach(LittleData)
            LittleData <- LittleData[which(Gene_name==Name),]
            detach(LittleData)
            #final<-results
        }
        else
        {
            stop(sprintf("Can't find gene: %s",Name))
        }  
    }

    ##add a little value to avoid 0 on denominators
    LittleData[,3]<-LittleData[c(0:nrow(LittleData)),3]+0.000001
    LittleData[,4]<-LittleData[c(0:nrow(LittleData)),4]+0.000001

    ##Calculate fold changes (both normal and logarithmyc values)
    LittleData$"FoldChanges"
    LittleData$"Log2FoldChanges"
    for(line in c(0:nrow(LittleData)))
    {
        FC=(LittleData[line,Tissue1]/LittleData[line,Tissue2])
        LittleData[line,"FoldChanges"]=FC       
        
        Log2FC=log2(FC)
        LittleData[line,"Log2FoldChanges"]=Log2FC       
    }
    
    return(LittleData)
}

Plotter<-function(Genes,Transcripts)
    ##from a vector of genes creates a table and a histogram
{
    Transcripts_FPKM<-Transcripts[,grep("FPKM",colnames(Transcripts))]
    Transcripts_FPKM$gene_name<-Transcripts$gene_name
    subsetted=NULL
    
    for (gene in Genes)
    {
        temp<-subset(Transcripts_FPKM,Transcripts_FPKM$gene_name==gene)
        subsetted<-rbind(subsetted,temp)	
    }

    subsetted_ok=NULL
    for (gene in Genes)
    {
        if (length(grep(gene,subsetted$gene_name)>1))
        {
            multiple<-subsetted[grep(gene,subsetted$gene_name,),]
            mean_vec<-apply(multiple[,-ncol(multiple)],2,mean)
            ##mean_vec$gene_nane
            ##mean_vec
            subsetted_ok<-rbind(subsetted_ok,mean_vec)
        }
        else
        {
            subsetted_ok<-rbind(subsetted_ok,subsetted[grep(gene,subsetted$gene_name,),-ncol(subsetted)])
        }
    }
    subsetted<-subsetted_ok
    ##subsetted$gene_name<-Genes
    row.names(subsetted)<-Genes
    colnames(subsetted)<-gsub("FPKM.","",colnames(subsetted))
    colnames(subsetted)<-gsub("\\."," ",colnames(subsetted))

    subsetted<-subsetted[,-ncol(subsetted)]
    subsetted<-as.data.frame(subsetted)

    tsubsetted<-t(subsetted)

    ##row.names(tsubsetted)<-colnames(subsetted)

    tsubsetted<-as.data.frame(tsubsetted)
    tsubsetted$tissue<-row.names(tsubsetted)
    
    for (gene in Genes)
    {
        out_file <- paste(gene,"png",sep=".")
        ylab_ok <- paste(gene,"(FPKM)",sep=" ")
        df2 <- tsubsetted[,c(gene, "tissue")]
        print(gene)
        p<-ggplot(data=df2)+
            geom_bar(stat="identity",aes_string(x="tissue",y=gene,fill="tissue"))+
            theme_minimal()+
            labs(x="TISSUE",y=ylab_ok)+
            theme(axis.text.x=element_text(angle=75,vjust=0.5,size=10))+
            theme(axis.title.x=element_text(size=15))+
            theme(axis.text.y=element_text(size=15))+
            theme(axis.title.y=element_text(size=15))+
            theme(legend.position="none")
        
        png(filename=out_file,width = 1200, height = 800)
        
        print(p)
        dev.off()
    }
}