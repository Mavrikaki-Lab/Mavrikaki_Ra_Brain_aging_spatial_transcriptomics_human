# 02_BayesSpace.R
# ------------------------------
# This script runs BayesSpace clustering on Visium data using Harmony-corrected PCA.
# It loops through different values of k (k=8 and k=16) and saves clustering results and visualizations.


# Load required libraries
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(spatialLIBD)
library(BayesSpace)
library(RColorBrewer)
library(ggplot2)
library(Polychrome)

# Define paths
work_dir <- here("02_BayesSpace")
plots_dir <- here("02_BayesSpace", "BS_Plots")
processed_data_dir <- here("01_BuildSpe", "BS_ProcessedData")

# Load filtered SPE object
load(file = file.path(processed_data_dir, "spe_filtered_final.Rdata"), verbose = TRUE)

# Preprocess SPE for BayesSpace (normalize + PCA)
spe <- spatialPreprocess(spe)
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# Define range of cluster numbers
Ks <- seq(8, 16)
set.seed(30122)

# Offset check visualization (to verify sample separation)
pdf(file = file.path(plots_dir, "BayesSpace_offset_check.pdf"))
clusterPlot(spe, "sample_id", color = NA) + 
  labs(fill = "Sample", title = "Offset check")
dev.off()

# Run BayesSpace clustering for each K using Harmony
for (k in Ks) {
  message("Running BayesSpace with k = ", k)
  spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)
  
  # Store clustering result with unique column name
  spe$bayesSpace_temp <- spe$spatial.cluster
  bayesSpace_name <- paste0("bayesSpace_harmony_", k)
  colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name
  
  # Save clustering results as CSV per sample
  cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(processed_data_dir, "clustering_results")
  )
  
  # Visualization: per-sample spatial clusters
  pdf(file = file.path(plots_dir, paste0("vis_clus_bayesSpace_harmony_", k, ".pdf")))
  for (i in seq_along(unique(spe$sample_id))) {
    my_plot <- vis_clus(
      spe = spe,
      clustervar = bayesSpace_name,
      sampleid = unique(spe$sample_id)[i],
      colors = setNames(Polychrome::palette36.colors(k), seq_len(k))
    )
    print(my_plot)
  }
  dev.off()
}

# Session info for reproducibility
session_info()
