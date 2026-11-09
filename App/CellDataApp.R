#First version created on 7/8/2025 by Vedant Desai 

#_____LIBRARIES_____#
library(tidyverse)
library(ggplot2)
library(tidyr)
library(shiny)
library(shinyWidgets)
library(bslib)
library(plotly)
library(bsicons)
library(qs2)
library(here)

# Load data files once globally to share across all sessions
HD_proc <- qs2::qs_read("HDproc.qs2")
V2_proc <- qs2::qs_read("V2proc.qs2")

#_____APP SECTION_____#

#Setting basic UI of the page
ui <- fluidPage(
  #formatting that helps when viewing on different platforms
  tags$style("
    body {
    -moz-transform: scale(0.9, 0.9); zoom: 0.9; zoom: 100%;
    }
  "),
  titlePanel(
    h1("Cell Type Analysis Plots", align="center")
  ), 
  page_navbar(
    #Page that displays the HD section
    nav_panel(
      title="HD",
      sidebarLayout(
        fluid=TRUE,
        sidebarPanel(
          width=3,
          selectizeInput("gene_input", "Gene Name", choice=NULL),
          checkboxGroupInput("plot_layers", "Select Plots:", c("Violin Plot"="violin", "Box Plot"="boxPlot", "Jitter Points"="jitterPoint"), selected="violin")
        ),
        mainPanel(
          width=9,
          card(
            card_body(plotOutput("cell_plot")),
            height=1000
          ),
          card(
            card_image(file=here("Acknowledgements.png")),
            card_body(
              hr(),
              p("Developed by Vedant Desai", align="center")
            )
          )
        )
      )
    ),
    #Page that displays the V2 section
    nav_panel(
      title="V2", 
      sidebarLayout(
        fluid=TRUE,
        sidebarPanel(
          width=3,
          selectizeInput("seq_input", "Seq", choice=NULL),
          checkboxGroupInput("plot_layers", "Select Plots:", c("Violin Plot"="violin", "Box Plot"="boxPlot", "Jitter Points"="jitterPoint"), selected="violin")
        ),
        mainPanel(
          width=9,
          card(
            card_body(plotOutput("V2_plot")),
            height=1000
          ),
          card(
            card_image(file=here("Acknowledgements.png")),
            card_body(
              hr(),
              p("Developed by Vedant Desai", align="center")
            )
          )
        )
      )
    )
  )
)

#______________________________________SERVER___________________________________

#Setting the server side of the app
server <- function(input, output, session){
  updateSelectizeInput(session, "gene_input", choice=unique(HD_proc$gene_name), server=TRUE, selected="VAMP1")
  updateSelectizeInput(session, "seq_input", choice=unique(V2_proc$Seq), server=TRUE)
  #filtering the original datasets to contain the cell and age specific information on the gene that the user inputted 
  filteredData <- reactive({
    req(input$gene_input)
    selected_gene <- input$gene_input
    
    cD_geneData <- HD_proc %>% 
      filter(gene_name==selected_gene)
    
    list(cell=cD_geneData)
  })
  
  V2filteredData <- reactive({
    req(input$seq_input)
    selected_gene <- input$seq_input
    
    V2_geneData <- V2_proc %>% 
      filter(Seq==selected_gene)
    
    list(seq=V2_geneData)
  })
  
  
  #Outputting the HD plot based on how the user wants to display their data 
  output$cell_plot <- renderPlot({
    cD_Data <- filteredData()$cell
    
    p <- ggplot(cD_Data, aes(x=cellType, y=Expression, fill=age)) +
      facet_wrap(~cellGroup, scales="free", nrow=3) +
      theme_minimal()+
      theme(
        panel.border = element_rect(color="black", fill=NA, linewidth=0.75),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1, face="bold"),
        axis.text.y=element_text(face="bold"),
        strip.text=element_text(face="bold")
      )
    
    if("violin" %in% input$plot_layers){
      p <- p + geom_violin(trim=FALSE, position=position_dodge(0.9), alpha=0.5)
    }
    if ("boxPlot" %in% input$plot_layers){
      p <- p + geom_boxplot(position=position_dodge(0.9), width=0.15, color="black")
    } 
    if("jitterPoint" %in% input$plot_layers){
      p <- p + geom_jitter(stat="identity", alpha=0.5, position=position_dodge(0.9))
    }
    p
  })
  
  #Outputting the V2 plot based on how the user wants to display their data 
  output$V2_plot <- renderPlot({
    V2_Data <- V2filteredData()$seq
    
    p2 <- ggplot(V2_Data, aes(x=Cell, y=Exp, fill=V2age)) +
      facet_wrap(~cellGroup, scales="free", nrow=3) +
      theme_minimal()+
      theme(
        panel.border = element_rect(color="black", fill=NA, linewidth=0.75),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1, face="bold"),
        axis.text.y=element_text(face="bold"),
        strip.text=element_text(face="bold")
      )
    
    if("violin" %in% input$plot_layers){
      p2 <- p2 + geom_violin(trim=FALSE, position=position_dodge(0.9), alpha=0.5)
    }
    if ("boxPlot" %in% input$plot_layers){
      p2 <- p2 + geom_boxplot(position=position_dodge(0.9), width=0.15, color="black")
    } 
    if("jitterPoint" %in% input$plot_layers){
      p2 <- p2 + geom_jitter(stat="identity", alpha=0.5, position=position_dodge(0.9))
    }
    p2
  })
}

#Displays the app
shinyApp(ui = ui, server = server)

