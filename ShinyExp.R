#First version created on 7/8/2025 by Vedant Desai 
#Purpose: Visualization of whole-genome expression data collected using Visium High Definition (2 um resolution) and Visium v2 (55 um resolution) from 
#   human dorsolateral prefrontal cortex (DLPFC) tissue of young (18-39 years old), middle (40-69 years old), and old (70-97 years old) age groups. 

#___________________________________LIBRARIES___________________________________

library(tidyverse)
library(readxl)
library(ggplot2)
library(tidyr)
library(shiny)
library(shinyWidgets)
library(bslib)
library(plotly)
library(bsicons)
library(here)

#___________________________________DATA PROCESSING SECTION_____________________

#Uploading Files and storing in list for processing 
V2 <- read_excel(here("pseudobulk_vst_by_domain_V2_Copy (3).xlsx"), sheet="Pool")
HD <- read_excel(here("pseudobulk_counts_vst_general_code_age_layer - Copy.xlsx"), sheet="in")
fileList <- list(HD, V2)

#Function that pivots the data from wide to long format for easier processing  
makeLongerFunc <- function(data){
  long_data <- data %>% 
    pivot_longer(
      cols=-gene_name,
      names_to="Sample",
      values_to="Expression"
    )
  return(long_data)
}
long_dataList <- lapply(fileList, makeLongerFunc)

#Extracting age and layer information from the long data of HD to be able to plot those characteristics in their specific groups 
HD_age <- "(?<=_)young|old|middle(?=-)"
HD_layer <- "(?<=-)$|L2|L1|L3|L4|L5|L6|WM"
HD_proc <- long_dataList[[1]] %>% 
  mutate(
    age=str_extract(Sample, HD_age),
    layer=str_extract(Sample, HD_layer)
  ) %>% 
  filter(!is.na(age)&!is.na(layer)) %>% 
  mutate(Expression=as.numeric(Expression)) %>% 
  mutate(
    layer=case_when(
      layer=="L3" ~ "L3_L4", 
      layer=="L4" ~ "L3_L4", 
      TRUE~layer
    )
  ) %>% 
  mutate(
    layer=factor(
      layer,
      levels=c("L1", "L2", "L3_L4", "L5", "L6", "WM")
    ),
    age=factor(age, levels=c("young", "middle", "old"))
  ) 

#Extracting age and layer information from the long data of V2 to be able to plot those characteristics in their specific groups 
V2_age <- "(?<=_)young|middle|old(?=_)"
V2_layer <- "Sp08D01|Sp08D02|Sp08D03|Sp08D04|Sp08D05|Sp08D06|Sp08D07|Sp08D08(?=_)"
V2_proc <- long_dataList[[2]] %>% 
  mutate(
    age=str_extract(Sample, V2_age),
    layer=str_extract(Sample, V2_layer),
    layer=recode(layer, Sp08D01="L1", Sp08D02="L5", Sp08D03="WM1", Sp08D04="L2/L3", Sp08D05="L3", Sp08D06="WM2", Sp08D07="L6", Sp08D08="L3/L4")
  ) %>% 
  mutate(Expression=as.numeric(Expression)) %>% 
  mutate(
    layer=factor(
      layer,
      levels=c("L1", "L2/L3", "L3", "L3/L4", "L5", "L6", "WM1", "WM2")
    ),
    age=factor(age, levels=c("young", "middle", "old"))
  ) 
procList <- list(HD_proc, V2_proc)

#Function that gets the median expression across age groups in each layer of both the HD and V2 datasets
medianExpFunction <- function(data){
  medianExp <- data %>% 
    group_by(age, layer, gene_name) %>% 
    summarize(
      medExp=median(Expression, na.rm=TRUE),
      .groups="drop"
    )
  plotData <- data %>% 
    full_join(medianExp, by=c("gene_name", "layer", "age"))
  return(plotData)
}
finDatalist <- lapply(procList, medianExpFunction)

#________________________________SHINY SECTION_______________________________________

#Setting the user interface of the app
ui <- fluidPage(
  titlePanel(
    h1("Aging-related Transcriptomic Changes with Spatial and Cell Type-Specific Resolution in the Human Brain", align="center")
  ), 
  card(
    card_header(
      h5("Purpose", align="center")
    ),
    card_body(
      p("This platform enables visualization of whole-genome expression data collected using spatial transcriptomics approach (Visium High-Definition (2 μm resolution) and
Visium v2 (55 μm resolution) from human dorsolateral prefrontal cortex (DLPFC) tissue of a young (18-39 years), middle (40-69 years), and old (70-97 years) group,
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
            accordion_panel("HD plot", plotOutput("HD_plot")), 
            accordion_panel("V2 plot", plotOutput("V2_plot"))
          )
        )
      ),
      card(
        card_header(
          h5("Acknowledgements", align="center")
        ), 
        card_image(file=here("Acknowledgements.png"))
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
      )
    )
  )
)

#Setting the server side of the app
server <- function(input, output, session){
  updateSelectizeInput(session, "multi_input", choices=unique(HD_proc$gene_name), server=TRUE)
  updateSelectizeInput(session, "gene_input", choice=unique(HD_proc$gene_name), server=TRUE, selected="VAMP1")
  
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
          facet_wrap(~layer, nrow=1) +
          labs(title=input$gene_input)+
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
          facet_wrap(~age, nrow=1)+
          labs(title=input$gene_input)+
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
      p
    })
    
    #Outputting the HD plot based on how the user wants to display their data 
    output$HD_plot <- renderPlot({
      HD_Data <- filteredData()$hd
      
      if("ageGroups" %in% input$x_axis){
        p2 <- ggplot(HD_Data, aes(x=age, y=Expression, fill=age, text=paste("Age:", age))) +
          facet_wrap(~layer, nrow=1) +
          labs(title=input$gene_input)+
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
          facet_wrap(~age, nrow=1)+
          labs(title=input$gene_input)+
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
      p2
    })
    
    #Both the HD and V2 heatmap outputs 
    output$HD_HM <- renderPlot({
      ggplot(data=heatMapFiltered()$hdHM, aes(x=layer, y=gene_name, fill=medExp))+
        geom_tile()+
        facet_wrap(~age) +
        scale_fill_gradient(low="lightblue", high="darkred", name="Median Expression")
    })
    output$V2_HM <- renderPlot({
      ggplot(data=heatMapFiltered()$v2HM, aes(x=layer, y=gene_name, fill=medExp))+
        geom_tile()+
        facet_wrap(~age)+
        scale_fill_gradient(low="lightblue", high="darkred", name="Median Expression")
    })
}

shinyApp(ui = ui, server = server)




