# 03_Pseudobulk.R
# ------------------------------
# This script generates pseudobulk data by BayesSpace domains (k = 8)
# and performs preliminary PCA + DGE modeling using spatialLIBD.


# Load libraries
library(here)
library(spatialLIBD)
library(BayesSpace)
library(RColorBrewer)
library(ggplot2)
library(sessioninfo)
library(scater)
library(jaffelab)

# Set BayesSpace clustering resolution
k <- 8

# Output directories
dir_rdata <- here("03_Pseudobulk", "PS_ProcessedData")
dir_plots <- here("03_Pseudobulk", "PS_Plots")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE object
spe_path <- here("01_BuildSpe", "BS_ProcessedData", "spe_filtered_final.Rdata")
stopifnot(file.exists(spe_path))
load(spe_path)

# Import BayesSpace clustering results
cluster_dir <- here("02_BayesSpace", "BS_ProcessedData", "clustering_results", paste0("bayesSpace_harmony_", k))
spe <- cluster_import(spe, cluster_dir = cluster_dir)

# Format BayesSpace domain variable
spe$BayesSpace <- factor(
  paste0("Sp", sprintf("%02d", k), "D", sprintf("%02d", spe[[paste0("imported_bayesSpace_harmony_", k)]]))
)
levels(spe$BayesSpace) <- c("Sp08D04", "Sp08D05", "Sp08D07", "Sp08D06", "Sp08D01", "Sp08D03", "Sp08D08", "Sp08D02")

# Run pseudobulk aggregation by domain and sample
sce_pseudo <- registration_pseudobulk(
  spe,
  var_registration = "BayesSpace",
  var_sample_id = "sample_id",
  min_ncells = 10
)
dim(sce_pseudo)

# Simplify metadata
colData(sce_pseudo) <- colData(sce_pseudo)[, sort(c(
  "sample_id", "BayesSpace", "sex", "tissue_bank", "group"
))]

# Save SCE pseudobulk
saveRDS(sce_pseudo, file = file.path(dir_rdata, paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")))

# Add variables for modeling functions
sce_pseudo$registration_variable <- sce_pseudo$BayesSpace
sce_pseudo$registration_sample_id <- sce_pseudo$sample_id

# Define modeling parameters
covars <- c("tissue_bank", "sex")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

# Fit registration model
registration_mod <- registration_model(sce_pseudo, covars = covars)
block_cor <- registration_block_cor(sce_pseudo, registration_model = registration_mod)

# Run differential expression tests
results_enrichment <- registration_stats_enrichment(sce_pseudo, block_cor, covars, gene_ensembl, gene_name)
results_pairwise <- registration_stats_pairwise(sce_pseudo, registration_mod, block_cor, gene_ensembl, gene_name)

# ANOVA test k ≥ 3
results_anova <- if (k >= 3) {
  registration_stats_anova(
    sce_pseudo, block_cor, covars, gene_ensembl, gene_name, suffix = suffix
  )
} else {
  NULL
}

# Save results
modeling_results <- list(
  "anova" = results_anova,
  "enrichment" = results_enrichment,
  "pairwise" = results_pairwise
)
save(modeling_results, file = file.path(dir_rdata, paste0("modeling_results_BayesSpace_k", sprintf("%02d", k), ".Rdata")))

# Reproducibility
Sys.time()
proc.time()
session_info()
