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
        if(TissueCol%in%colnames(Transcripts))
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

SearchByTissue<-function(Tissue,Name,Transcripts)
    ##subsets the dataframe to extract only columns reguarding a certain tissue
{
    print(Tissue)
    Name<-toupper(Name)
    print(Name)
    ##create FPKM.tissue and OUT.tissue strings
    paste("Plotting FKPM for tissue ",Tissue)
    Fpkm="FPKM"
    Out="OUT"

    TissueCol=paste(Fpkm,Tissue,"sep"=".")
    ##valutate the value of TissueCol and search it in Transcripts
    temp<-quote(TissueCol)

    #same with cov
    Cov="cov"
    TissueCov=paste(Cov,Tissue,"sep"=".")
    temp<-quote(TissueCov)

    TranscriptsInfos<-Transcripts[2:10]

    final<-Transcripts[2:10]
    final$FPKM<-Transcripts[,TissueCol]
    final$cov<-Transcripts[,TissueCov]

    ##rename the columns
    colnames(final)<-c("Chromosome","Strand","Start","End","t_name","Num_exons","Length","Gene_ID","Gene_Name","FPKM","Coverage")
    
    if (!is.null(Name) & !is.na(Name) & Name!="")
    {

        if(Name %in% final$Gene_Name)
        {
            attach(final)
            final <- final[which(Gene_Name==Name),]
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
    print("Search By Gene")
    Name<-toupper(Name)
    print(Name)

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
    print("Search By Gene Feature")
    Name<-toupper(Name)
    print(Name)
    print(Feature)
    
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
    print("Search By Differential Expression")
    print(Tissue1)
    print(Tissue2)
    Name<-toupper(Name)
    print(Name)
    
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
    
    LittleData<-Transcripts[,c(6,10)]
    LittleData[,3]<-TransCov1
    LittleData[,4]<-TransCov2

    colnames(LittleData)<-c("ID","Gene_Name",Tissue1,Tissue2)

    if (!is.null(Name) && !is.na(Name) && Name!="")
    {
        print(!is.null(Name))
        print(!is.na(Name))
        print(Name)
    }
    
    if (!is.na(Name) & !is.null(Name) & Name!="")
    {
        if(Name %in% LittleData$Gene_Name)
        {
            attach(LittleData)
            LittleData <- LittleData[which(Gene_Name==Name),]
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
                                        #lapply(Genes,toupper())
    Genes<-toupper(Genes)
    print(Genes)
    
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

    if (length(Genes)==1)
    {
        print("subsetted is ok")
    }
    else
    {
        subsetted<-subsetted[,-ncol(subsetted)]
        subsetted<-as.data.frame(subsetted)
    }
    tsubsetted<-t(subsetted)

    row.names(tsubsetted)<-colnames(subsetted)

    tsubsetted<-as.data.frame(tsubsetted)
    tsubsetted$tissue<-row.names(tsubsetted)
    
    for (gene in Genes)
    {
        out_file <- paste(gene,"png",sep=".")
        ylab_ok <- paste(gene,"(FPKM)",sep=" ")
        
        #print(tsubsetted)
        if (length(Genes)==1)
        {
            df2<-tsubsetted
        }
        else
        {
            df2 <- tsubsetted[,c(gene, "tissue")]
        }
        print(gene)
        p<-ggplot(data=df2)+
            geom_bar(stat="identity",aes_string(x="tissue",y=gene,fill="tissue"))+
            theme_minimal()+
            labs(x="TISSUE",y=ylab_ok)+
            theme(axis.text.x=element_text(angle=75,vjust=0.5,size=15))+
            theme(axis.title.x=element_text(size=15))+
            theme(axis.text.y=element_text(size=15))+
            theme(axis.title.y=element_text(size=15))+
            theme(legend.position="none")
        
    }
    return(p)
}

SearchTranscriptGroup<-function(query_type,Name,ModulesFile,Transcripts)
{
    ##Name: gene name to be analyzed (es gene77) or gene symbol (es STK33)
    ##ModulesFile: file containig informations about coexpression
    ##modules (TOM files from WGCNA pipeline)

    load(ModulesFile)
    
    if(query_type=="symbol")
    {
        #all gene symbols are upper cases
        Name<-toupper(Name)
        
        colnames(datExpr)<-as.character(Transcripts[match(colnames(datExpr),transcripts$gene_id),]$gene_name)
        pos<-match(Name,colnames(datExpr))
    }

    ##probably redundant
    if(query_type=="gene_name")
    {
        pos<-match(Name,colnames(datExpr))
    }
    
    
    geneModule<-moduleLabels[pos]

    ##extracting all genes belongin to the module
    Gene_Module<-names(datExpr)[moduleLabels==geneModule]

    if(query_type=="symbol")
    {
        subsetted<-transcripts[transcripts$gene_name %in% Gene_Module,]
    }

    if(query_type=="gene_name")
    {
        subsetted<-transcripts[transcripts$gene_id %in% Gene_Module,]
    }
            
    return(subsetted[,c(1:10)])
}

Network<-function(query_type,query,MODULES,Transcripts,corr,results)
{
    ##file hash for output
    File_Hash<-paste(format(Sys.time(), "%Y%m%d%H%M"),"png",sep=".")
    
    ## == Sanity checks ==
    if(query_type!="gene_name" & query_type!="symbol")
    {
        print("Choose between 'gene_name' and 'symbol'")
        stop("Argument 1 not recognized")
    }

    if(file.exists(MODULES))
    {
        load(MODULES)
    }else{
        print("File not found. Makes sure that the file containing the modules is correct")
    }
    
    ##== END ==

    ## == Default vaules of correlation and number of results ==
    if(length(corr)==0)
    {
        corr=0.9
    }

    if(length(results)==0)
    {
        results=50
    }
    ## == End of Defaults ==

    ## importing transcripts data
    
    ##== The query is the gene symbol or gene name? ==
    
    if(query_type=="symbol")
    {
        #all gene symbols are upper cases
        query<-toupper(query)
        
        colnames(datExpr)<-as.character(transcripts[match(colnames(datExpr),transcripts$gene_id),]$gene_name)

        query_position<-match(query,colnames(datExpr))
        query_module<-moduleLabels[query_position]
    }

    ## probably redundant
    if(query_type=="gene_name")
    {
        query_position<-match(query,colnames(datExpr))
        query_module<-moduleLabels[query_position]
    }

    ##== END ==

    module_line<-which(moduleLabels %in% query_module)
    datExpr_module<-datExpr[,module_line]

    A<-cor(datExpr_module)
    print("Correlation matrix created")
    colnames(A)

    diag(A)<-0
    A[abs(A)<corr]=0

    row_col=match(query,colnames(A))
    cor_vec=abs(A[,row_col])
    C<-names(sort(abs(A[,row_col]),decreasing=TRUE)[c(0:results)])

    if(query%in%C){
        print("Query gene is in C")
    }else{
        C[results]<-query
    }
    D<-A[C,C]

    ##D/40 to stress distances between nodes
    D<-D/40

    print("Create graph")

    
    if(query_type=="symbol")
    {
        graph<-graph_from_adjacency_matrix(D,mode="max",weighted=TRUE,diag=FALSE)

        graph<-simplify(graph)

        
        V(graph)$color<-"white"
        V(graph)[query]$color<-"red"

        out_file=paste("Network",query,File_Hash,sep="_")
        ##plotting
        png(out_file,width=1080,heigh=720)
        plot(graph)
                                        #plot(graph,vertex.shape="vrectangle")
        dev.off()
    }
    
    if(query_type=="gene_name")
    {

        graph<-graph_from_adjacency_matrix(D,mode="max",weighted=TRUE,diag=FALSE)
            
        V(graph)$gene=as.character(transcripts[match(C,transcripts$gene_id),]$gene_name)
        #V(graph)$gene=as.character(transcripts[match(C,transcripts$t_name),]$gene_name)
        coul=rainbow(n=length(V(graph)$gene))
        my_color=coul[as.numeric(as.factor(V(graph)$gene))]
        
        V(graph)$color<-"white"
        V(graph)[query]$color<-"red"
        V(graph)$shape<-"vrectangle"
        graph<-simplify(graph)

        ##plotting

        out_file=paste("Network",query,File_Hash,sep="_")
        
        png(out_file,width=1080,heigh=720)
        
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
        dev.off()
    }    
    return(graph)
}
