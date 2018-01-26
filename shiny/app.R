library(shiny)

#Before running this chunk of code, import and format correctly the data



#!!!! Attention: keep the "(" of fluid page on the same line
ui <- fluidPage(
    titlePanel("Search By Tissue"),
    
    sidebarLayout
    (
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
        mainPanel
        (
            tableOutput(outputId="SBT")
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
}

shinyApp(ui, server)
