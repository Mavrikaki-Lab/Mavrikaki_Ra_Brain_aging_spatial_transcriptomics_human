#First version created on 7/8/2025 by Vedant Desai 
#Purpose: Visualization of whole-genome expression data collected using Visium High Definition (2 um resolution) and Visium v2 (55 um resolution) from 
#   human dorsolateral prefrontal cortex (DLPFC) tissue of young (18-39 years old), middle (40-69 years old), and old (70-97 years old) age groups. 

#___________________________________LIBRARIES___________________________________

library(tidyverse)
library(ggplot2)
library(tidyr)
library(shiny)
library(shinyWidgets)
library(bslib)
library(plotly)
library(bsicons)
# library(ggprism)
# library(grafify)
library(qs2)
library(here)

#___________________________________GLOBAL DATA LOADING SECTION_____________________
# Load data files once globally to share across all sessions
finDatalist <- qs2::qs_read("finDatalist.qs2")
HD_proc <- qs2::qs_read("HD_proc.qs2")
V2_proc <- qs2::qs_read("V2_proc.qs2")

# Define regex patterns for age and layer extraction
V2_age <- "(?<=_)young|middle|old(?=_)"
V2_layer <- "Sp08D01|Sp08D02|Sp08D03|Sp08D04|Sp08D05|Sp08D06|Sp08D07|Sp08D08(?=_)"
HD_age <- "(?<=_)young|old|middle(?=-)"
HD_layer <- "(?<=-)$|L2|L1|L3|L4|L5|L6|WM"

#________________________________SHINY SECTION_______________________________________

#Setting the user interface of the app
ui <- fluidPage(
  title = "Brain Aging Spatial Transcriptomics",
  tags$style("
    body {
    -moz-transform: scale(0.9, 0.9); zoom: 0.9; zoom: 100%;
    }
  "),
  titlePanel(
    h1("Aging-related Transcriptomic Changes with Spatial and Cell Type-Specific Resolution in the Human Brain", align="center")
  ), 
  card(
    card_header(
      h5("Overview", align="center")
    ),
    card_body(
      p("This platform enables visualization of whole-genome expression data collected using spatial transcriptomics approach (Visium High-Definition (16 μm resolution) and
Visium V2 (55 μm resolution)) from human dorsolateral prefrontal cortex (DLPFC) tissue of a young (18-39 years), middle (40-69 years), and old (70-97 years) group,
allowing characterization of age-related molecular changes within and across cortical layers (and white matter), and cell types. ​
From Rai et al., to be submitted", align="center")
    )
  ),
  page_navbar(
    nav_panel(
      title="Individual Gene Expression Plots",
      sidebarLayout(
        sidebarPanel(
          width=3,
          selectizeInput("gene_input", "Gene Name", choice=NULL), 
          checkboxGroupInput("plot_layers", "Select Plots:", c("Violin Plot"="violin", "Box Plot"="boxPlot", "Jitter Points"="jitterPoint"), selected="violin"),
          radioGroupButtons("x_axis", "Select X axis:", c("Age Group"="ageGroups", "Layers"="brainLayers"), selected="brainLayers")
        ),
        mainPanel(
          width=9,
          accordion(
            open=c("HD plot", "V2 plot"),
            accordion_panel(
              "HD plot", 
              accordion(
                height=300,
                accordion_panel("HD cluster", div(imageOutput("HDcluster"), style="text-align:center;"))
              ),
              plotOutput("HD_plot")
            ), 
            accordion_panel(
              "V2 plot", 
              accordion(
                height=275,
                accordion_panel("V2 cluster", div(imageOutput("V2cluster"), style="text-align:center;"))
              ),
              plotOutput("V2_plot")
            )
          )
        )
      ),
      card(
        card_image(file=here("Acknowledgements.png")),
        card_body(
          hr(),
          p("Developed by Vedant Desai", align="center")
        )
      )
    ),
    nav_panel(
      title="Multiple Gene Heatmaps",
      navset_card_underline(
        id="dSet",
        title="Datasets",
        sidebar=sidebar(
          h4("Select Genes"),
          selectizeInput("multi_input", "Genes to Graph", choices=NULL, multiple=TRUE)
        ),
        nav_panel("HD", plotOutput("HD_HM")), 
        nav_panel("V2", plotOutput("V2_HM"))
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

#Setting the server side of the app
server <- function(input, output, session){
  updateSelectizeInput(session, "multi_input", choices=unique(HD_proc$gene_name), server=TRUE)
  updateSelectizeInput(session, "gene_input", choice=unique(HD_proc$gene_name), server=TRUE, selected="VAMP1")
 
  age_colors <- c(young = "#D22B2B", middle = "#0000FF", old = "#4F7942")
  
  output$HDcluster <- renderImage({
    list(src=here("HDclusters.png"), width="90%", height="225")
    },
    deleteFile=FALSE 
  )
  
  output$V2cluster <- renderImage({
    list(src=here("V2clusters.png"), width="75%", height="200")
   },
   deleteFile=FALSE
  )
  
  #filtering the original datasets to contain the layer and age specific information on the gene that the user inputted 
  filteredData <- reactive({
    req(input$gene_input)

    selected_gene <- input$gene_input
    
    HD_geneData <- finDatalist[[1]] %>% 
      filter(gene_name==selected_gene)
    V2_geneData <- finDatalist[[2]] %>% 
      filter(gene_name==selected_gene)
    
    list(hd=HD_geneData, v2=V2_geneData)
  })
  
  #filtering the original datasets to contain the layer and age specific information on the select genes the user inputted 
  heatMapFiltered <- reactive({
    req(input$multi_input)
    selected_gene <- input$multi_input
    
    HM_HDgeneData <- finDatalist[[1]] %>% 
      subset(gene_name %in% input$multi_input)
    HM_V2geneData <- finDatalist[[2]] %>% 
      subset(gene_name %in% input$multi_input)
    
    list(hdHM=HM_HDgeneData, v2HM=HM_V2geneData)
  })
    
    #Outputting the V2 plot based on how the user wants to display their data 
  output$V2_plot <- renderPlot({
    V2_Data <- filteredData()$v2
      
    if("ageGroups" %in% input$x_axis){
      p <- ggplot(V2_Data, aes(x=age, y=Expression, fill=age, group=age, text=paste("Age:", age))) +
        facet_grid(.~layer) +
        xlab("Age Groups")
      if("violin" %in% input$plot_layers) {
        p <- p + geom_violin(trim=FALSE, alpha=0.5)
      }
      if("boxPlot" %in% input$plot_layers) {
        p <- p + geom_boxplot(width=0.15, color="black")
      }
      if("jitterPoint" %in% input$plot_layers){
        p <- p + geom_jitter(stat="identity", position=position_jitterdodge(jitter.width=0.2), alpha=0.5)
      }
    }
    if("brainLayers" %in% input$x_axis){
      p <- ggplot(V2_Data, aes(x=layer, y=Expression, fill=age, text=paste("Layer:", layer))) +
        facet_grid(.~age, labeller=labeller(age=c(young="Young", middle="Middle", old="Old"))) +
        xlab("Layers")
      if("violin" %in% input$plot_layers) {
        p <- p + geom_violin(trim=FALSE, alpha=0.5)
      }
      if("boxPlot" %in% input$plot_layers) {
        p <- p + geom_boxplot(width=0.15, color="black")
      }
      if("jitterPoint" %in% input$plot_layers){
        p <- p + geom_jitter(stat="identity", position=position_jitterdodge(jitter.width=0.2), alpha=0.5)
      }
    }
    p <- p + labs(title=input$gene_input) +
      scale_fill_manual(values=age_colors) +
      theme_minimal()+
      theme(
        panel.border=element_rect(color="black", fill=NA, linewidth=0.75),
        strip.text=element_text(face="bold")
      )
    p
  })
    
    #Outputting the HD plot based on how the user wants to display their data 
output$HD_plot <- renderPlot({
    HD_Data <- filteredData()$hd
    if("ageGroups" %in% input$x_axis){
      p2 <- ggplot(HD_Data, aes(x=age, y=Expression, fill=age, text=paste("Age:", age))) +
        facet_grid(.~layer) +
        xlab("Age Groups") 
      if("violin" %in% input$plot_layers) {
        p2 <- p2 + geom_violin(trim=FALSE, alpha=0.5)
      }
      if("boxPlot" %in% input$plot_layers) {
        p2 <- p2 + geom_boxplot(width=0.15, color="black")
      }
      if("jitterPoint" %in% input$plot_layers){
        p2 <- p2 + geom_jitter(stat="identity", position=position_jitterdodge(jitter.width=0.2), alpha=0.5)
      }
    }
    if("brainLayers" %in% input$x_axis){
      p2 <- ggplot(HD_Data, aes(x=layer, y=Expression, fill=age, text=paste("layer:", layer))) +
        facet_grid(.~age, labeller=labeller(age=c(young="Young", middle="Middle", old="Old")))+
        xlab("Layers") 
      if("violin" %in% input$plot_layers) {
        p2 <- p2 + geom_violin(trim=FALSE, alpha=0.5)
      }
      if("boxPlot" %in% input$plot_layers) {
        p2 <- p2 + geom_boxplot(width=0.15, color="black")
      }
      if("jitterPoint" %in% input$plot_layers){
        p2 <- p2 + geom_jitter(stat="identity", position=position_jitterdodge(jitter.width=0.2), alpha=0.5)
      }
    }
    p2 <- p2 + labs(title=input$gene_input) +
      scale_fill_manual(values=age_colors) +
      theme_minimal()+
      theme(
        panel.border=element_rect(color="black", fill=NA, linewidth=0.75),
        strip.text=element_text(face="bold")
      ) 
    p2
  })
    
    #Both the HD and V2 heatmap outputs 
    output$HD_HM <- renderPlot({
      ggplot(data=heatMapFiltered()$hdHM, aes(x=layer, y=gene_name, fill=medExp))+
        geom_tile()+
        facet_grid(.~age, labeller=labeller(age=c(young="Young", middle="Middle", old="Old"))) +
        scale_fill_gradientn(colors=c("blue", "white", "red"), name="Median Expression") +
        theme_minimal()+
        theme(
          panel.border=element_rect(color="black", fill=NA, linewidth=0.75),
          strip.text=element_text(face="bold")
        ) +
        xlab("Layers") + ylab("Gene Name")
        
    })
    output$V2_HM <- renderPlot({
      ggplot(data=heatMapFiltered()$v2HM, aes(x=layer, y=gene_name, fill=medExp))+
        geom_tile()+
        facet_grid(.~age, labeller=labeller(age=c(young="Young", middle="Middle", old="Old")))+
        scale_fill_gradientn(colors=c("blue", "white", "red"), name="Median Expression")+
        theme_minimal()+
        theme(
          panel.border=element_rect(color="black", fill=NA, linewidth=0.75),
          strip.text=element_text(face="bold")
        ) +
        xlab("Layers") + ylab("Gene Name")
    })
}

shinyApp(ui = ui, server = server)




