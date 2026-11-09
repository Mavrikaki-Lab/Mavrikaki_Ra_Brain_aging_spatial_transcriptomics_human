# 01_BuildSpe.R
# ------------------------------
# This script builds a SpatialExperiment (SPE) object from 10X Visium output.
# Includes metadata integration, mitochondrial content calculation,
# filtering of spots and genes, QC with scran, HVG selection, PCA, UMAP, t-SNE,
# batch correction with Harmony, and preparation for BayesSpace.


# Load required libraries
library(here)
library(readxl)
library(sessioninfo)

library(SpatialExperiment)
library(spatialLIBD)
library(scuttle)
library(scran)
library(scater)        
library(PCAtools)      
library(BiocParallel)
library(harmony)

# Define paths
work_dir <- here("01_BuildSpe")
plots_dir <- here("01_BuildSpe", "BS_Plots")
processed_data_dir <- here("01_BuildSpe", "BS_ProcessedData")
metadata_file <- here("01_BuildSpe", "BS_Scripts", "MetaDataV2_Anupama.xlsx")

# Read metadata
df_metadata <- read_excel(metadata_file, sheet = "MasterSheet")

# Build sample information table
sample_info <- data.frame(
  sample_id = df_metadata$SampleID,
  sample = df_metadata$Sample,
  group = df_metadata$Group,
  group_label = df_metadata$Group_label,
  age = df_metadata$Age,
  age_label = df_metadata$Age_label,
  age_bin = df_metadata$AgeBin,
  sex = df_metadata$Sex,
  tissue_bank = df_metadata$TissueBank,
  race = df_metadata$Race,
  ethnicity = df_metadata$Ethnicity,
  PMI_hours = df_metadata$PMI_hours,
  manner_death = df_metadata$MannerDeath,
  date_sectioned = df_metadata$DateSectioned,
  date_processed = df_metadata$DateProcessed,
  days_post_sectioning = df_metadata$DaysPostSectioning,
  processing_batch = df_metadata$ProcessingBatch,
  SI_cycle = df_metadata$SIcycle,
  bank_RIN = df_metadata$BankRIN,
  measured_RIN = df_metadata$MeasuredRIN,
  neuropathology_diagnosis = df_metadata$NeuropathalogyDiagnosis
)

# Add and validate sample paths
sample_info$sample_path <- file.path(here("01_BuildSpe", "SR_samples"), sample_info$sample_id, "outs")
stopifnot(all(file.exists(sample_info$sample_path)))
write.csv(sample_info, file = file.path(processed_data_dir, "spe_combined.csv"), row.names = FALSE)

# Build SPE object
cat("Processing and combining the samples...\n")
spe <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = here("reference", "refdata-gex-GRCh38-2020-A", "genes", "genes.gtf")
)

# Add sample metadata to colData
colData(spe)$sample_id <- spe$sample_id
metadata_vars <- setdiff(colnames(sample_info), c("sample_path"))
for (var in metadata_vars) {
  colData(spe)[[var]] <- sample_info[[var]][match(spe$sample_id, sample_info$sample_id)]
}

# Add segmentation info (tissue_positions.csv)
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
  file <- file.path(sample_info$sample_path[sample_info$sample_id == sampleid], "spatial", "tissue_positions.csv")
  if (!file.exists(file)) return(NULL)
  x <- read.csv(file, header = FALSE)
  colnames(x) <- c("barcode", "in_tissue", "row", "col", "imagerow", "imagecol")
  x$key <- paste0(x$barcode, "_", sampleid)
  return(x)
})
segmentations <- do.call(rbind, segmentations_list)
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <- segmentations[segmentation_match, ]
segmentation_info <- segmentation_info[, !colnames(segmentation_info) %in% c("barcode", "row", "col", "imagerow", "imagecol", "key")]
colData(spe) <- cbind(colData(spe), segmentation_info)

# Mitochondrial expression
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

# Save full object before filtering
save(spe, file = file.path(processed_data_dir, "spe_before_qc.RData"))

# Filter genes with no expression
no_expr <- which(rowSums(counts(spe)) == 0)
spe <- spe[-no_expr, ]

# Keep only spots inside tissue
spe <- spe[, colData(spe)$in_tissue == 1]

# Visualization of QC metrics
vis_grid_clus(spe_raw, clustervar = "in_tissue",
              pdf = file.path(plots_dir, "in_tissue_grid.pdf"),
              sort_clust = FALSE,
              colors = c("TRUE" = "grey90", "FALSE" = "orange"))

# Scran QC
qcstats <- perCellQCMetrics(spe, subsets = list(Mito = is_mito))
qcfilter <- quickPerCellQC(qcstats, sub.fields = "subsets_Mito_percent")

# Store flags in colData
spe$scran_discard <- factor(qcfilter$discard)
spe$scran_low_lib_size <- factor(isOutlier(spe$sum_umi, type = "lower", log = TRUE, batch = spe$sample_id))
spe$scran_low_n_features <- factor(qcfilter$low_n_features)
spe$scran_high_subsets_Mito_percent <- factor(qcfilter$high_subsets_Mito_percent)

# Save plots
for (i in colnames(qcfilter)) {
  vis_grid_clus(
    spe = spe,
    clustervar = paste0("scran_", i),
    pdf = file.path(plots_dir, paste0("scran_", i, ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
  )
}

# Remove low-library spots
spe <- spe[, colData(spe)$scran_low_lib_size == FALSE]

# Save final filtered object
save(spe, file = file.path(processed_data_dir, "spe_filtered_final.RData"))

# Normalize with scran
set.seed(123)
spe$scran_quick_cluster <- quickCluster(spe, BPPARAM = MulticoreParam(3), block = spe$sample_id)
spe <- computeSumFactors(spe, clusters = spe$scran_quick_cluster, BPPARAM = MulticoreParam(3))
spe <- logNormCounts(spe)

# Model variance
dec <- modelGeneVar(spe, block = spe$sample_id, BPPARAM = MulticoreParam(3))
top.hvgs <- getTopHVGs(dec, prop = 0.1)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
save(top.hvgs, top.hvgs.fdr5, top.hvgs.fdr1, file = file.path(processed_data_dir, "top.hvgs_all.RData"))

# PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)

# Elbow plot
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
pdf(file.path(plots_dir, "pca_elbow_final.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

# UMAP + TSNE
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
spe <- runTSNE(spe, dimred = "PCA", name = "TSNE_perplexity80", perplexity = 80)

# Harmony
spe <- RunHarmony(spe, "sample_id", verbose = FALSE)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")

# BayesSpace adjustments
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- colData(spe)$array_col
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# Final save
save(spe, file = file.path(processed_data_dir, "spe_filtered_final.RData"))

# Reproducibility info
Sys.time()
proc.time()
session_info()
