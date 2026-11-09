library(rstudioapi)
library(tidyverse)
setwd(fs::path_dir(getSourceEditorContext()$path))
# NOTE: This code will check version, this is our recommendation, it may work
#.      other versions
stopifnot(R.version$major>= 4) # requires R4
# stopifnot(compareVersion(R.version$minor,"3.3")==0) # requires >=4.3.3
stopifnot(compareVersion(as.character(BiocManager::version()), "3.18")>=0)
stopifnot(compareVersion(as.character(packageVersion("Seurat")), "5.0.0")>=0)

.libPaths('/n/app/bcbio/R4.3.1_singlecell/')
library(Seurat)
library(data.table)
library(qs)

### Set up run information
data_dir <- "../data/"

samples <- list.files(data_dir)
samples <- samples[samples != "1756" & samples != '1512_D8']

################################################################################
# make seurat obj for 8 um raw 
################################################################################

### Make individual seurat objects for each sample
for (i in 1:length(samples)){
  print(i)
  seurat_obj <- Load10X_Spatial(paste0(data_dir, samples[i]), 
                                filename = 'raw_feature_bc_matrix.h5',
                                slice = samples[i],
                                bin.size = c(8))
  assign(paste0(samples[i], "_seurat"),
         seurat_obj) # stores Seurat object in variable of corresponding sample name
}

### Merge all seurat objects

seurat_ID <- paste0(samples, "_seurat") # get names of all objects


u <- get(seurat_ID[2])
for (i in 3:length(seurat_ID)) {
  u <- c(u, get(seurat_ID[i]))
} ## makes a list of all seurat objects

seurat_merge <- merge(x = get(seurat_ID[1]),
                      y = u,
                      add.cell.ids = samples,
                      project = "visium_hd")

qsave(seurat_merge, '../data_processed/seurat_merged.qs')

seurat_merge <- qread('../data_processed/seurat_merged.qs')

# Mitochondrial genes for mouse genome
idx <- grep("^MT-", rownames(GetAssay(seurat_merge, "Spatial.008um")))
rownames(GetAssay(seurat_merge, "Spatial.008um"))[idx][1:10]
# Mitochondrial genes vs. nuclear genes ratio
seurat_merge$mitoRatio <- PercentageFeatureSet(object = seurat_merge, pattern = "^MT-")
seurat_merge$mitoRatio <- seurat_merge@meta.data$mitoRatio/100 # Divide by 100 for Ratio instead of Percentage

# Number of genes per UMI for each cell
seurat_merge$Log10GenesPerUMI <- log10(seurat_merge$nFeature_Spatial.008um) / log10(seurat_merge$nCount_Spatial.008um)
seurat_merge$barcode <- rownames(seurat_merge@meta.data)

# Extract cell level metadata
metadata <- seurat_merge@meta.data
# metadata$barcode <- rownames(metadata)

# # Import experimental metadata
# metaexp <- read.csv("/path/to/experimental/metadata/meta.csv")
# 
# # Check matching of IDs
# all(metaexp$sample %in% metadata$orig.ident)
# all(metadata$orig.ident %in% metaexp$sample)
# 
# #change headings to match
# colnames(metaexp)[1] <- "orig.ident"
# 
# metafull <- plyr::join(metadata, metaexp,
#                        by = c("orig.ident"))
# 
# # Replace seurat object metadata
# if(all(metafull$barcode == rownames(seurat_merge@meta.data))) {
#   rownames(metafull) <- metafull$barcode
#   seurat_merge@meta.data <- metafull
# }


## Join layers (each sample is a separate layer)
seurat_merge[["Spatial.008um"]] <- JoinLayers(seurat_merge[["Spatial.008um"]])

### Save Seurat object for future processing
qsave(seurat_merge, file = "../data_processed/seurat_pre-QC.qs")
# write.csv(seurat_merge@meta.data, file = "metadata_pre-filtered.csv")

seurat_merge <- qread('../data_processed/seurat_pre-QC.qs')

################################################################################
# make seurat obj for 8 um spaceranger-filtered
################################################################################

### Make individual seurat objects for each sample
for (i in 1:length(samples)){
  print(i)
  seurat_obj <- Load10X_Spatial(paste0(data_dir, samples[i]), 
                                slice = samples[i],
                                bin.size = c(8))
  assign(paste0(samples[i], "_seurat"),
         seurat_obj) # stores Seurat object in variable of corresponding sample name
}

### Merge all seurat objects

seurat_ID <- paste0(samples, "_seurat") # get names of all objects


u <- get(seurat_ID[2])
for (i in 3:length(seurat_ID)) {
  u <- c(u, get(seurat_ID[i]))
} ## makes a list of all seurat objects

seurat_merge <- merge(x = get(seurat_ID[1]),
                      y = u,
                      add.cell.ids = samples,
                      project = "visium_hd")

qsave(seurat_merge, '../data_processed/seurat_filtered_merged.qs')

seurat_merge <- qread('../data_processed/seurat_filtered_merged.qs')

# Mitochondrial genes for mouse genome
idx <- grep("^MT-", rownames(GetAssay(seurat_merge, "Spatial.008um")))
rownames(GetAssay(seurat_merge, "Spatial.008um"))[idx][1:10]
# Mitochondrial genes vs. nuclear genes ratio
seurat_merge$mitoRatio <- PercentageFeatureSet(object = seurat_merge, pattern = "^MT-")
seurat_merge$mitoRatio <- seurat_merge@meta.data$mitoRatio/100 # Divide by 100 for Ratio instead of Percentage

# Number of genes per UMI for each cell
seurat_merge$Log10GenesPerUMI <- log10(seurat_merge$nFeature_Spatial.008um) / log10(seurat_merge$nCount_Spatial.008um)
seurat_merge$barcode <- rownames(seurat_merge@meta.data)

# Extract cell level metadata
metadata <- seurat_merge@meta.data
# metadata$barcode <- rownames(metadata)

# # Import experimental metadata
# metaexp <- read.csv("/path/to/experimental/metadata/meta.csv")
# 
# # Check matching of IDs
# all(metaexp$sample %in% metadata$orig.ident)
# all(metadata$orig.ident %in% metaexp$sample)
# 
# #change headings to match
# colnames(metaexp)[1] <- "orig.ident"
# 
# metafull <- plyr::join(metadata, metaexp,
#                        by = c("orig.ident"))
# 
# # Replace seurat object metadata
# if(all(metafull$barcode == rownames(seurat_merge@meta.data))) {
#   rownames(metafull) <- metafull$barcode
#   seurat_merge@meta.data <- metafull
# }


## Join layers (each sample is a separate layer)
seurat_merge[["Spatial.008um"]] <- JoinLayers(seurat_merge[["Spatial.008um"]])

### Save Seurat object for future processing
qsave(seurat_merge, file = "../data_processed/seurat_filtered_pre-QC.qs")
# write.csv(seurat_merge@meta.data, file = "metadata_pre-filtered.csv")

seurat_merge <- qread('../data_processed/seurat_filtered_pre-QC.qs')

################################################################################
# make seurat obj for 16 um raw 
################################################################################

### Make individual seurat objects for each sample
for (i in 1:length(samples)){
  print(i)
  seurat_obj <- Load10X_Spatial(paste0(data_dir, samples[i]), 
                                filename = 'raw_feature_bc_matrix.h5',
                                slice = samples[i],
                                bin.size = c(16))
  seurat_obj$orig.ident <- samples[i]
  assign(paste0(samples[i], "_seurat"),
         seurat_obj) # stores Seurat object in variable of corresponding sample name
}

### Merge all seurat objects

seurat_ID <- paste0(samples, "_seurat") # get names of all objects


u <- get(seurat_ID[2])
for (i in 3:length(seurat_ID)) {
  u <- c(u, get(seurat_ID[i]))
} ## makes a list of all seurat objects

seurat_merge <- merge(x = get(seurat_ID[1]),
                      y = u,
                      add.cell.ids = samples,
                      project = "visium_hd")

qsave(seurat_merge, '../data_processed/seurat_16um_merged.qs')

seurat_merge <- qread('../data_processed/seurat_16um_merged.qs')

# Mitochondrial genes for mouse genome
idx <- grep("^MT-", rownames(GetAssay(seurat_merge, "Spatial.016um")))
rownames(GetAssay(seurat_merge, "Spatial.016um"))[idx][1:10]
# Mitochondrial genes vs. nuclear genes ratio
seurat_merge$mitoRatio <- PercentageFeatureSet(object = seurat_merge, pattern = "^MT-")
seurat_merge$mitoRatio <- seurat_merge@meta.data$mitoRatio/100 # Divide by 100 for Ratio instead of Percentage

# Number of genes per UMI for each cell
seurat_merge$Log10GenesPerUMI <- log10(seurat_merge$nFeature_Spatial.016um) / log10(seurat_merge$nCount_Spatial.016um)
seurat_merge$barcode <- rownames(seurat_merge@meta.data)

# Extract cell level metadata
metadata <- seurat_merge@meta.data
# metadata$barcode <- rownames(metadata)

# # Import experimental metadata
# metaexp <- read.csv("/path/to/experimental/metadata/meta.csv")
# 
# # Check matching of IDs
# all(metaexp$sample %in% metadata$orig.ident)
# all(metadata$orig.ident %in% metaexp$sample)
# 
# #change headings to match
# colnames(metaexp)[1] <- "orig.ident"
# 
# metafull <- plyr::join(metadata, metaexp,
#                        by = c("orig.ident"))
# 
# # Replace seurat object metadata
# if(all(metafull$barcode == rownames(seurat_merge@meta.data))) {
#   rownames(metafull) <- metafull$barcode
#   seurat_merge@meta.data <- metafull
# }


## Join layers (each sample is a separate layer)
seurat_merge[["Spatial.016um"]] <- JoinLayers(seurat_merge[["Spatial.016um"]])

### Save Seurat object for future processing
qsave(seurat_merge, file = "../data_processed/seurat_16um_pre-QC.qs")
# write.csv(seurat_merge@meta.data, file = "metadata_pre-filtered.csv")

seurat_merge <- qread('../data_processed/seurat_16um_pre-QC.qs')
################################################################################

ImageFeaturePlot(seurat_merge, 'nCount_Spatial.008um', max.cutoff = 'q80')
ImageFeaturePlot(seurat_merge, 'nCount_Spatial.008um', max.cutoff = 'q80', fov = 'Hct15HAAA038_D5.008um')
