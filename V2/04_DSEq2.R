# 04_DESeq2.R
# ------------------------------
# This script performs domain-specific differential expression using DESeq2.
# It compares "Old" and "Young" vs "Middle" across 8 BayesSpace domains (k = 8),
# applying LFC shrinkage with ashr and exporting Excel results per comparison.
# Author: Irais Valenzuela

# Load libraries
library(here)
library(SingleCellExperiment)
library(SpatialExperiment)
library(DESeq2)
library(tidyverse)
library(ashr)
library(writexl)
library(stringr)

# Define paths
input_dir <- here("04_DESeq2", "DS_Input")
output_dir <- here("04_DESeq2", "DS_Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load pseudobulk SCE and gene mapping
sce <- readRDS(file.path(input_dir, "sce_pseudo_BayesSpace_k08_age_layer.rds"))
mapping <- readRDS(file.path(input_dir, "ens_symbol_mapping.rds"))  # gene_id, gene_name

# Create DESeq2 object
dds_full <- DESeqDataSetFromMatrix(
  countData = counts(sce),
  colData = colData(sce),
  design = ~ age_layer + sex
)

# Filter genes with low counts
dds_full <- dds_full[rowSums(counts(dds_full)) > 10, ]
dds_full <- DESeq(dds_full)

# Define domains and comparisons
layers <- paste0("Sp08D0", 1:8)
comparisons <- list(
  c("Old", "Middle"),
  c("Young", "Middle")
)

# Run DE analysis with LFC shrinkage
de_results_shrunk <- list()

for (layer in layers) {
  for (comparison in comparisons) {
    group1 <- paste0(comparison[1], "_", layer)
    group2 <- paste0(comparison[2], "_", layer)
    
    if (all(c(group1, group2) %in% levels(dds_full$age_layer))) {
      message("Running: ", group1, " vs ", group2)
      
      res_shrunk <- lfcShrink(dds_full,
                              contrast = c("age_layer", group1, group2),
                              type = "ashr")
      
      # Format results
      res_df <- as.data.frame(res_shrunk) %>%
        rownames_to_column("gene_id") %>%
        left_join(mapping, by = "gene_id") %>%
        mutate(
          svalue = ifelse("svalue" %in% colnames(.), svalue, NA),
          layer = layer,
          comparison = paste(group1, "vs", group2)
        ) %>%
        relocate(gene_id, gene_name, svalue)
      
      result_name <- paste0(layer, "_", comparison[1], "_vs_", comparison[2])
      de_results_shrunk[[result_name]] <- res_df
    } else {
      message("Skipping ", group1, " vs ", group2, ": missing samples")
    }
  }
}

# Organize results by comparison
split_comparisons <- function(results_list) {
  comparison_map <- list("Old_vs_Middle" = list(), "Young_vs_Middle" = list())
  for (name in names(results_list)) {
    comp <- results_list[[name]]$comparison[1]
    if (str_detect(comp, "^Old_.*vs Middle$")) {
      comparison_map$Old_vs_Middle[[name]] <- results_list[[name]]
    } else if (str_detect(comp, "^Young_.*vs Middle$")) {
      comparison_map$Young_vs_Middle[[name]] <- results_list[[name]]
    }
  }
  return(comparison_map)
}

shrunk_results_by_comparison <- split_comparisons(de_results_shrunk)

# Function to create summary tables
make_summary <- function(result_list) {
  lapply(result_list, function(df) {
    data.frame(
      padj_0.05 = sum(!is.na(df$padj) & df$padj < 0.05),
      strong_deg = sum(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1),
      comparison = unique(df$comparison),
      layer = unique(df$layer)
    )
  }) %>% bind_rows()
}

# Export Excel per comparison
for (comparison in names(shrunk_results_by_comparison)) {
  result_list <- shrunk_results_by_comparison[[comparison]]
  if (length(result_list) > 0) {
    summary_df <- make_summary(result_list)
    
    excel_data <- c(
      lapply(result_list, function(x) {
        sheetname <- paste0("Shrunk_", unique(x$layer))
        setNames(list(x), sheetname)
      }) %>% unlist(recursive = FALSE),
      list(Summary = summary_df)
    )
    
    write_xlsx(excel_data, path = file.path(output_dir, paste0("DEGs_Shrunk_", comparison, ".xlsx")))
  }
}

# Optional: save RDS
saveRDS(de_results_shrunk, file = file.path(output_dir, "deseq2_shrunk_results_filtered.rds"))

# Reproducibility info
cat("Reproducibility information:\n")
Sys.time()
sessionInfo()
