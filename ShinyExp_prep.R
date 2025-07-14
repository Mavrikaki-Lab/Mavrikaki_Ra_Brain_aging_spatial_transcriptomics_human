library(tidyverse)
library(qs2)

#Uploading Files and storing in list for processing 
V2 <- read_csv(here("pseudobulk_vst_by_domain_V2.csv"))
HD <- read_csv(here("pseudobulk_counts_vst_general_code_age_layer.csv"))
fileList <- list(HD, V2)

# #Function that pivots the data from wide to long format for easier processing  
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
qs2::qs_save(long_dataList, "long_dataList.qs2")

#Extracting age and layer information from the long data of HD to be able to plot those characteristics in their specific groups 
HD_age <- "(?<=_)young|old|middle(?=-)"
HD_layer <- "(?<=-)$|L2|L1|L3|L4|L5|L6|WM"
HD_proc <- long_dataList[[1]] %>%
  mutate(
    age=str_extract(Sample, HD_age),
    layer=str_extract(Sample, HD_layer), 
    layer=recode(layer, L2="L2/L3")
  ) %>%
  filter(!is.na(age)&!is.na(layer)) %>%
  mutate(Expression=as.numeric(Expression)) %>%
  mutate(
    layer=case_when(
      layer=="L3" ~ "L3/L4",
      layer=="L4" ~ "L3/L4",
      TRUE~layer
    )
  ) %>%
  mutate(
    layer=factor(
      layer,
      levels=c("L1", "L2/L3", "L3/L4", "L5", "L6", "WM")
    ),
    age=factor(age, levels=c("young", "middle", "old"))
  )
qs2::qs_save(HD_proc, "HD_proc.qs2")


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
  ) %>% 
  filter(!is.na(age)&!is.na(layer))

procList <- list(HD_proc, V2_proc)
qs2::qs_save(V2_proc, "V2_proc.qs2")

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
qs2::qs_save(finDatalist, "finDatalist.qs2")
