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
            sidebarLayout(
                sidebarPanel(
                    ## sliding bar that sets a value
                    textInput
                    (
                        inputId="gene",
                        label="Gene you want to analyze",
                        value="STK33"
                    ),
                    actionButton
                    (
                        inputId="search_plot",
                        label="Search"
                    )
                ),
                mainPanel(
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
                        inputId="tissue",
                        label="Tissue Of Interest",
                        choices=phenodata$tissue,
                        selected="rumen"
                    ),
                    textInput
                    (
                        inputId="gene",
                        label="Gene to be searched",
                        value="STK33",
                        placeholder="Your gene"
                    ),
                    actionButton
                    (
                        inputId="search",
                        label="Search"
                    )
                ),
                mainPanel(
                    tableOutput(outputId="SBT")
                )
            )
        )
    )
)

## Define server logic to plot various variables against mpg ----
server <- function(input, output) {
    output$SBT<-renderTable({
        input$search
        
        TISSUE<-isolate(input$tissue)
        GENE<-isolate(input$gene)
        SearchByTissue(TISSUE,GENE,transcripts)
    })
    
    output$plot<-renderPlot({
        input$search_plot
        GENE<-isolate(input$gene)
        Plotter(GENE,transcripts)
    })
}

shinyApp(ui, server)
