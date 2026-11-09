#First version created on 1/4/2026 by Vedant Desai 

#_____LIBRARIES_____#
library(tidyverse)
library(here)
library(qs2)

#_____DATA PROCESSING_____#

#Uploading HD and V2 datasets
HD <- read_csv(here("pseudobulk_counts_vst_general_code_age_cell_type.csv"))
V2 <- read_csv(here("cell2location_q05_cell_abundance_by_barcode.csv"))

#Pivoting HD data to long format for easier processing
HD_long <- HD %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "Sample", 
    values_to = "Expression"
  )

#Pivoting V2 data to long format for easier processing
V2_long <- V2 %>% 
  pivot_longer(
    cols = -key,
    names_to = "Cell",
    values_to = "Exp"
  )

#Variables that processing step looks for when extracting relevant information
V2_age <- "(?<=_)$|V2AG[0-9]+"
V2_young <- "(?<=_)$|V2AG2|V2AG37|V2AG3|V2AG1|V2AG10|V2AG12|V2AG4|V2AG5|V2AG6|V2AG7|V2AG8|V2AG9"
V2_cellGroupPattern <- "Exc|In(?=_)"

#Processing step: extracts age and cell group information from the V2 dataset
V2_proc <- V2_long %>% 
  mutate(V2age=str_extract(key, V2_age)) %>%
  mutate(cellGroup=str_extract(Cell, V2_cellGroupPattern),
         cellGroup=replace_na(cellGroup, "Glial")) %>% 
  mutate(Seq=str_remove(key, V2_age)) %>% 
  mutate(V2age = case_match(V2age, c("V2AG2", "V2AG37", "V2AG3","V2AG1", "V2AG10", "V2AG12", "V2AG4", "V2AG5", "V2AG6", "V2AG7", "V2AG8", "V2AG9") ~ "young",
                            c("V2AG13", "V2AG14", "V2AG15", "V2AG16", "V2AG17", "V2AG18", "V2AG19", "V2AG20", "V2AG21", "V2AG22", "V2AG23") ~ "middle", 
                            .default="old")) %>% 
  mutate(
    V2age=factor(
      V2age, levels=c("young", "middle", "old")
    )
  ) %>% 
  filter(!is.na(V2age)&!is.na(cellGroup)&!is.na(Exp))

#Variables that processing step looks for when extracting relevant information
HD_age <- "(?<=_)young|middle|old(?=-)"
HD_layer <- "(?<=-)Astro-FB|Astro-PP|Exc-L2-3|Exc-L3-5|Exc-L4-6-1|Exc-L4-6-2|Exc-L4-6-3|Exc-L5-6-1|Exc-L5-6-2|Exc-L5-6-HTR2C|In-LAMP5-1|In-LAMP5-2|In-PVALB-Ba|In-PVALB-Ch|In-RELN|In-SST|In-VIP|Microglia|Oligodendrocyte|OPC|Unknown"
HD_cellGroupPattern <- "(?<=-)$|Exc|In(?=-)"

#Processing step: extracts age, cell type, and cell group information from the HD dataset
HD_proc <- HD_long %>% 
  mutate(age=str_extract(Sample, HD_age), 
         cellType=str_extract(Sample, HD_layer),
         cellGroup=str_extract(Sample, HD_cellGroupPattern),
         cellGroup=replace_na(cellGroup, "Glial")) %>% 
  mutate(Expression=as.numeric(Expression)) %>% 
  filter(!is.na(age)&!is.na(cellType)&!is.na(Expression)) %>% 
  mutate(
    cellType=factor(
      cellType,
      levels=c("Astro-FB", "Astro-PP", "Exc-L2-3", "Exc-L3-5", "Exc-L4-6-1", "Exc-L4-6-2", "Exc-L4-6-3", "Exc-L5-6-1", "Exc-L5-6-2", "Exc-L5-6-HTR2C", "In-LAMP5-1", "In-LAMP5-2", "In-PVALB-Ba", "In-PVALB-Ch", "In-RELN", "In-SST", "In-VIP", "Microglia", "Oligodendrocyte", "OPC", "Unknown")
    ),
    age=factor(age, levels=c("young", "middle", "old"))
  ) 

#Saves processed information into files for app
qs2::qs_save(HD_proc, "HDproc.qs2")
qs2::qs_save(V2_proc, "V2proc.qs2")

