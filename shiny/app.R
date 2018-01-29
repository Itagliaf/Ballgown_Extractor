library(shiny)
library(ggplot2)
##=== Loading functions ====
source("definitions.R")

##==== Loading and formatting data ====
## load("data.fil.RData")
## transcripts<-data.fil@expr$trans
## head(transcripts)
## phenodata<-read.csv("phenodata.csv","header"=TRUE)
## colnames(transcripts)<-NameFormatter(transcripts,phenodata)
##==== END ====

#!!!! Attention: keep the "(" of fluid page on the same line
ui <- fluidPage( #Begins a flui page (change when windows dimensions change)
    navlistPanel( #Starts a panel environment: each search is a panel
        tabPanel(
            #A better name must be found
            "Plotter",
            verticalLayout
            (
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
                ),
                plotOutput(outputId="plot")
               
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
            sidebarLayout(
                sidebarPanel(
                    ## sliding bar that sets a value
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
                mainPanel(
                    tableOutput(outputId="SBG")
                )
            )
        ),
        tabPanel(#Definition of a tab called search by gene feture
            "Search By Gene Feature",
            sidebarLayout(
                sidebarPanel(
                    ## sliding bar that sets a value
                    textInput
                    (
                        inputId="gene_feature",
                        label="Gene to be searched",
                        value="STK33",
                        placeholder="Your gene"
                    ),
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
                mainPanel(
                    tableOutput(outputId="SBGF")
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
    
    ## #output search by tissue

    output$SBT<-renderTable({
        input$search_tissue
        TISSUE<-isolate(input$tissue_tissue)
        GENE<-isolate(input$gene_tissue)
        SearchByTissue(TISSUE,GENE,transcripts)
    })

    output$SBG<-renderTable({
        input$search_gene
        
        GENE<-isolate(input$gene_gene)
        SearchByGene(GENE,transcripts)
    })
    
    ##output search by gene feature
    output$SBGF<-renderTable({
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
}

shinyApp(ui, server)
