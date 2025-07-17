# SpatialRegistration.R
# ------------------------------
# This script correlates BayesSpace domain-level enrichment stats
# with cortical layer modeling results using spatialLIBD.
# It uses top 100 genes by t-statistics and generates a correlation matrix.

# Load libraries
library(spatialLIBD)
library(tidyverse)
library(xlsx)
library(jaffelab)
library(ComplexHeatmap)
library(here)
library(circlize)
library(sessioninfo)

# Load custom plotting function
source(here("05_SpatialRegistration", "SR_Scripts", "layer_stat_cor_plot_custom.R"))

# Define directories
input_dir <- here("SpatialRegistration", "SR_ProcessedData")
plot_dir <- here("SpatialRegistration", "SR_Plots")

# Load registration results (currently only for k = 8)
k_list <- c(8)

# Load modeling results per k
registration_files <- map(k_list, ~ here(
  input_dir,
  paste0("annotated_modeling_results_BayesSpace_k", sprintf("%02d", .x), ".Rdata")
))
registration_data <- lapply(registration_files, function(x) get(load(x)))

# Extract enrichment t-statistics
registration_t_stats <- map(registration_data, function(data) {
  x <- data$enrichment
  t_stats <- x[, grep("^t_stat_", colnames(x))]
  colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))  # Remove "t_stat_" prefix
  return(t_stats)
})

# Quick check
map(registration_t_stats, jaffelab::corner)

# Load reference layer modeling results
layer_modeling_results <- fetch_data(type = "modeling_results")

# Compute correlation with layer model
cor_top100 <- map(
  registration_t_stats,
  ~ layer_stat_cor(
    .x,
    layer_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
  )
)

# Save correlation matrices
save(cor_top100, file = here(input_dir, "annotated_bayesSpace_layer_cor_top100.Rdata"))

# Optional: reorder for consistent plotting
cor_top100[[1]] <- cor_top100[[1]][order(rownames(cor_top100[[1]])), ]

# Plot correlation heatmaps using custom palette
pdf(here(plot_dir, "cor_top100_spatial_registration_customPalette.pdf"))
map(cor_top100, ~ layer_stat_cor_plot_custom(.x))
dev.off()

# Reproducibility info
Sys.time()
session_info()
