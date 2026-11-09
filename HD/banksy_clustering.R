# .libPaths('/n/app/bcbio/R4.4.0_spatial_dev/')
library(tidyverse)
library(patchwork)
library(Seurat)
library(qs)
library(bcbioR)
library(SeuratWrappers)
# library(SpatialExperiment)
library(Banksy)
options(future.globals.maxSize= 5000000000)

# this script needs ~150GB of memory to run

filtered <- qread('data_processed/seurat_16um_post-QC.qs')

################## BANKSY, one sample ##########################################

# see if banksy clustering removes some of the noise of seurat sketch clustering
# answer: yes, not perfectly but better

filtered_sample <- subset(filtered, subset = orig.ident == '2721_D7_Manual')

filtered_sample <- RunBanksy(filtered_sample, lambda = 0.2, verbose = T,
                             assay = 'Spatial.016um', slot = 'data', k_geom = 15)
filtered_sample <- RunPCA(filtered_sample, assay = "BANKSY", 
                          reduction.name = "pca.banksy", 
                          features = rownames(filtered_sample), npcs = 30)
filtered_sample <- FindNeighbors(filtered_sample, reduction = "pca.banksy", 
                                 dims = 1:30)
filtered_sample <- FindClusters(filtered_sample, cluster.name = "banksy_cluster",
                                resolution = 0.65)

color_pal = Seurat::DiscretePalette(n = length(unique(filtered_sample$banksy_cluster)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(filtered_sample$banksy_cluster))

SpatialDimPlot(filtered_sample, group.by = "banksy_cluster", 
               pt.size.factor = 7, cols = color_pal)+
  guides(fill=guide_legend(ncol=2, override.aes = list(size=3, alpha=1)))

#################### banksy for 2 samples, SPE #################################

# try to run banksy on multiple samples as the Banksy docs recommend (using SpatialExperiment objects)
# results: considerably less noisy than Seurat but took ~3 hours and 170GB memory 
# to run on 2 samples, not scaleable to all 16

samples <- c('2721_D7_Manual', '3192Y_D2')
images <- c('X2721_D7_Manual.016um', 'X3192Y_D2.016um')
seu_list <- lapply(samples, function(sample_name) {subset(filtered, subset = orig.ident == sample_name)})
seu_list <- lapply(seu_list, function(seu) {subset(seu, cells = rownames(GetTissueCoordinates(seu, scale = "hires")))})

seu_list <- lapply(seu_list, function(x) {
  NormalizeData(x, scale.factor = 5000, normalization.method = 'RC')
})

#' Compute HVGs
hvgs <- lapply(seu_list, function(x) {
  VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)

seurat_to_spe <- function(seu, sample_id, img_id) {
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(GetTissueCoordinates(seu, scale = "hires")[, c('x', 'y')])
  
  seu$sample_id <- NULL
  
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}

spe_list <- lapply(1:length(seu_list), function(i){seurat_to_spe(seu_list[[i]], samples[i], images[i])})

aname <- "normcounts"
spe_list <- Map(function(spe, seu) {
  assay(spe, aname) <- GetAssayData(seu)
  spe[hvgs,]
}, spe_list, seu_list)


compute_agf <- FALSE
k_geom <- 8
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, 
                   compute_agf = compute_agf, k_geom = k_geom)

spe_joint <- do.call(cbind, spe_list)

lambda <- 0.2
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)

res <- 0.65
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)

seu_joint <- as.Seurat(spe_joint)

qsave(seu_joint, 'data_processed/seurat_16um_two_samples_banksy.qs')

seu_joint@images <- c(seu_list[[1]]@images, seu_list[[2]]@images)

color_pal = Seurat::DiscretePalette(n = length(unique(seu_joint$clust_M0_lam0.2_k50_res0.65)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(seu_joint$clust_M0_lam0.2_k50_res0.65))
SpatialDimPlot(seu_joint, group.by = 'clust_M0_lam0.2_k50_res0.65', cols = color_pal)


####################### banksy for 2 samples, seurat ###########################

# run banksy on multiple samples as described in SeuratWrappers documentation
# results: not as clean as SPE Banksy, but still a little cleaner than sketch clustering

samples <- c('2721_D7_Manual', '3192Y_D2')
# samples <- unique(filtered$orig.ident)[1:6]
filtered_n2 <- subset(filtered, subset = orig.ident %in% samples)

coords <- lapply(names(filtered_n2@images), function(image) GetTissueCoordinates(filtered_n2, image, scale = "hires")) %>% bind_rows()
filtered_n2 <- AddMetaData(filtered_n2, coords)

filtered_n2 <- NormalizeData(filtered_n2)
filtered_n2 <- FindVariableFeatures(filtered_n2)

meta <- filtered_n2@meta.data %>% filter(!is.na(x), !is.na(y))
filtered_n2 <- subset(filtered_n2, cells = rownames(meta))


filtered_n2 <- RunBanksy(filtered_n2, assay = 'Spatial.016um', lambda = 0.8, 
                         k_geom = 8, group = 'orig.ident', slot = 'data', 
                         dimx = 'x', dimy = 'y', split.scale = TRUE)

filtered_n2 <- RunPCA(filtered_n2, assay = "BANKSY", 
                      reduction.name = "pca.banksy", 
                      features = rownames(filtered_n2), npcs = 30)
filtered_n2 <- FindNeighbors(filtered_n2, reduction = "pca.banksy", 
                             dims = 1:30)
filtered_n2 <- FindClusters(filtered_n2, cluster.name = "banksy_cluster",
                            resolution = 0.3)


color_pal = Seurat::DiscretePalette(n = length(unique(filtered_n2$banksy_cluster)),
                                    palette = "polychrome")

names(color_pal) <- sort(as.numeric(unique(filtered_n2$banksy_cluster)))

SpatialDimPlot(filtered_n2, group.by = "banksy_cluster", pt.size.factor = 7,
                                        cols = color_pal,
               ncol = 3) & 
  guides(fill = guide_legend(override.aes = list(size=3, alpha=1), ncol = 2))

qsave(filtered_n2, 'data_processed/seurat_16um_two_samples_banksy_seurat.qs')

####################### banksy for all samples, seurat #########################

# try to run banksy using SeuratWrappers instructions on all samples. need to use
# fewer genes (750) to avoid a vector being longer than R allows

coords <- lapply(names(filtered@images), function(image) GetTissueCoordinates(filtered, image, scale = "hires")) %>% bind_rows()
filtered <- AddMetaData(filtered, coords)

filtered <- NormalizeData(filtered)
filtered <- FindVariableFeatures(filtered, nfeatures = 750)

meta <- filtered@meta.data %>% filter(!is.na(x), !is.na(y))
filtered_subset <- subset(filtered, cells = rownames(meta))

filtered_subset <- subset(filtered_subset, features = VariableFeatures(filtered_subset))


filtered_subset <- RunBanksy(filtered_subset, assay = 'Spatial.016um', lambda = 0.8, 
                         k_geom = 8, group = 'orig.ident', slot = 'data', 
                         dimx = 'x', dimy = 'y', split.scale = TRUE, 
                         features = 'all')

# filtered_subset <- SketchData(
#   object = filtered_subset,
#   assay = 'BANKSY',
#   ncells = 50000,
#   method = "LeverageScore",
#   sketched.assay = "sketch_banksy",
#   features = VariableFeatures(filtered_subset)
# )
# 
# DefaultAssay(filtered_subset) <- "BANKSY"


gc()
filtered_subset <- RunPCA(filtered_subset, assay = "BANKSY", 
                      reduction.name = "pca.banksy", 
                      features = rownames(filtered_subset), npcs = 30)
filtered_subset <- FindNeighbors(filtered_subset, reduction = "pca.banksy", 
                             dims = 1:30)
filtered_subset <- FindClusters(filtered_subset, cluster.name = "banksy_cluster_res0.2",
                            resolution = 0.2)
gc()
filtered_subset <- FindClusters(filtered_subset, cluster.name = "banksy_cluster_res0.3",
                                resolution = 0.3)
filtered_subset <- FindClusters(filtered_subset, cluster.name = "banksy_cluster_res0.1",
                                resolution = 0.1)

# color_pal = Seurat::DiscretePalette(n = length(unique(filtered_subset$banksy_cluster_res0.1)),
#                                     palette = "polychrome")
# names(color_pal) <- sort(as.numeric(unique(filtered_subset$banksy_cluster_res0.1)))
# 
# SpatialDimPlot(filtered_subset, group.by = "banksy_cluster_res0.1", pt.size.factor = 7,
#                cols = color_pal, ncol = 6)&
#   guides(fill = guide_legend(override.aes = list(size=3, alpha=1), ncol = 2))

qsave(filtered_subset, 'data_processed/seurat_16um_banksy_seurat.qs')

