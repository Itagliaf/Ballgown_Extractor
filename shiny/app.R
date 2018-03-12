#!!!! Attention: keep the "(" of fluid page on the same line
ui <- fluidPage( #Begins a flui page (change when windows dimensions change)
    tabsetPanel(
        tabPanel(
            ##A better name must be found
            "Plotter",
           
            fluidRow
            (
                column
                (
                    width=4,
                    textInput
                    (
                        inputId="gene_plot",
                        label="Gene you want to analyze",
                        value="STK33"
                    ),
                    actionButton
                    (
                        inputId="search_plot",
                        label="Search"
                    )
                )
            ),
            fluidRow
            (
                column
                (
                    width=12,
                    plotOutput(outputId="plot")
                )
            )
        ),
        tabPanel( #Definition of a tab called search by tissue
            "Search By Tissue",
            sidebarLayout(
                sidebarPanel(
                    ## sliding bar that sets a value
                    selectInput
                    (
                        inputId="tissue_tissue",
                        label="Tissue Of Interest",
                        choices=phenodata$tissue,
                        selected="rumen"
                    ),
                    textInput
                    (
                        inputId="gene_tissue",
                        label="Gene to be searched",
                        placeholder="Your gene",
                        value="STK33"
                    ),
                    actionButton
                    (
                        inputId="search_tissue",
                        label="Search"
                    )
                ),
                mainPanel(
                    tableOutput(outputId="SBT")
                )
            )
        ),
        tabPanel(#Definition of a tab called search by gene
            "Search By Gene",
            fluidRow
            (
                column
                (
                    width=12,
                    textInput
                    (
                        inputId="gene_gene",
                        label="Gene to be searched",
                        placeholder="Your gene",
                        value="STK33"
                    ),
                    actionButton
                    (
                        inputId="search_gene",
                        label="Search"
                    )
                ),
                column
                (
                    width=12,
                    dataTableOutput(outputId="SBG")
                )
            )
        ),
        tabPanel(#Definition of a tab called search by gene feture
            "Search By Gene Feature",
            fluidRow(
                column(
                    width=12,
                    textInput
                    (
                        inputId="gene_feature",
                        label="Gene to be searched",
                        value="STK33",
                        placeholder="Your gene"
                    ),
                    #buttons to choose the feature (intron or exon)
                    radioButtons
                    (
                        inputId="feature_feature",
                        choices=c("exon","intron"),
                        selected="exon",
                        label="Feutered to be analyzed"
                    ),
                    actionButton
                    (
                        inputId="search_feature",
                        label="Search"
                    )
                ),
                column
                (
                    width=12,
                    dataTableOutput(outputId="SBGF")
                )
            )
        ),
        tabPanel(#Definition of a tab called search by differential fold
            "Search By Differential Fold",
            sidebarLayout(
                sidebarPanel(
                    ## sliding bar that sets a value
                    selectInput
                    (
                        inputId="tissue1_fold",
                        label="First Tissue to be analyzed",
                        choices=phenodata$tissue,
                        selected="liver"
                    ),
                    selectInput
                    (
                        inputId="tissue2_fold",
                        label="Second Tissue to be analyzed",
                        choices=phenodata$tissue,
                        selected="tongue"
                    ),
                    textInput
                    (
                        inputId="gene_fold",
                        label="Gene you want to analyze",
                        value="STK33"
                    ),
                    actionButton
                    (
                        inputId="search_fold",
                        label="Search"
                    )
                ),
                mainPanel(
                    tableOutput(outputId="SBDF")
                )
            )
        ),
        tabPanel(
            "Search Transcriptiona Module",
            sidebarLayout(
                sidebarPanel(
                    textInput
                    (
                        inputId="gene_module",
                        label="Gene you want to analyze",
                        value="gene1"
                    ),
                    actionButton
                    (
                        inputId="search_module",
                        label="Search"
                    )
                    
                ),
                mainPanel(
                    tableOutput(outputId="MO")
                )
            )
        ),
        tabPanel(
            "Plot Network Graph",           
            fluidRow
            (
                column
                (
                    width=4,
		    selectInput
                    (
                        inputId="gene_id_network",
                        label="Gene name or gene symbol available?",
                        choices=c("symbol","gene_name"),
                        selected="symbol"
                    ),
                    textInput
                    (
                        inputId="transcript_network",
                        label="Transcript you want to analyze",
                        value="STK3"
                    ),
                    selectInput
                    (
                        inputId="correlation_network",
                        label="Correlation threshold",
                        choices=c(0.75,0.8,0.85,0.9,0.95),
                        selected=0.8
                    ),
                    textInput
                    (
                        input="number_network",
                        label="Number of nodes to be plotted",
                        value=25
                    ),
                    actionButton
                    (
                        inputId="crate_network",
                        label="Create"
                    )
                )
            ),
            fluidRow
            (
                column
                (
                    width=12,
                    plotOutput
                    (
                        outputId="NT",
                        width="1080px",
                        height="720px",
			
                    )
                )
            )
        ),
	tabPanel(#Definition of a tab called search by gene feture
            "Gene fold in all tissues",
            fluidRow(
                column(
                    width=12,
                    textInput
                    (
                        inputId="gene_fold_all",
                        label="Gene to be searched",
                        value="STK33",
                        placeholder="Your gene"
                    ),
                    actionButton
                    (
                        inputId="search_fold_all",
                        label="Search"
                    )
                ),
                column
                (
                    width=12,
                    tableOutput(outputId="GFA")
                )
            )
        )
    )
)

## Define server logic to plot various variables against mpg ----
server <- function(input, output) {
    ##output plotter
    output$plot<-renderPlot({
        input$search_plot
        GENE<-isolate(input$gene_plot)
        Plotter(GENE,transcripts)
    })
    
    ##output search by tissue
    output$SBT<-renderTable({
        input$search_tissue
        TISSUE<-isolate(input$tissue_tissue)
        GENE<-isolate(input$gene_tissue)
        SearchByTissue(TISSUE,GENE,transcripts)
    })

    output$SBG<-renderDataTable({
        input$search_gene
        
        GENE<-isolate(input$gene_gene)
        SearchByGene(GENE,transcripts)
    })
    
    ##output search by gene feature
    output$SBGF<-renderDataTable({
        input$search_feature
        
        GENE<-isolate(input$gene_feature)
        FEATURE<-isolate(input$feature_feature)
        SearchByFeature(GENE,FEATURE,data.fil,phenodata)
    })
    
    output$SBDF<-renderTable({
        input$search_fold
        
        TISSUE1<-isolate(input$tissue1_fold)
        TISSUE2<-isolate(input$tissue2_fold)
        GENE<-isolate(input$gene_fold)
        SearchByDiffFoldExpr(TISSUE1,TISSUE2,GENE,transcripts)
    })

    output$MO<-renderTable({
        input$search_module
        GENE<-isolate(input$gene_module)
        SearchTranscriptGroup(GENE,"OK_Gene_Modules.Rdata","Little_1000_TOM-block.1.RData")
    })

    output$NT<-renderPlot({
        input$crate_network
	ID<-isolate(input$gene_id_network)
	TRAN<-isolate(input$transcript_network)
        COR<-isolate(input$correlation_network)
        NUM<-isolate(input$number_network)
	
        Network(ID,TRAN,"OK_Gene_Modules.Rdata",data.fil,COR,NUM)        
    })

    output$GFA<-renderTable({
        input$gene_fold_all
        GENE_FOLD<-isolate(input$gene_fold_all)
        GeneFoldTissue(GENE_FOLD,transcripts)
    })

}

shinyApp(ui, server)
