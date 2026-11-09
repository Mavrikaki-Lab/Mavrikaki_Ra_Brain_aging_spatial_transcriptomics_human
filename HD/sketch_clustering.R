.libPaths('/n/app/bcbio/R4.4.0_spatial_dev/')
library(tidyverse)
library(patchwork)
library(Seurat)
library(qs)
library(bcbioR)
library(SeuratWrappers)
library(SpatialExperiment)
library(Banksy)
options(future.globals.maxSize= 5000000000)

# this script needs ~150GB of memory to run

filtered <- qread('data_processed/seurat_16um_post-QC.qs')

################# no normalization, all samples ################################

# concern that normalization removes relevant biological signal from less 
# transcriptionally active cell types, so attempt clustering without normalization
# result: clustering noisy, driven by low info bins

filtered <- SetAssayData(filtered, assay = 'Spatial.016um', layer = 'data', 
                         new.data = filtered@assays$Spatial.016um@layers$counts)

filtered <- FindVariableFeatures(filtered)
filtered <- ScaleData(filtered)

filtered <- SketchData(
  object = filtered,
  assay = 'Spatial.016um',
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

qsave(filtered, 'data_processed/seurat_16um_sketched.qs')

filtered <- qread('data_processed/seurat_16um_sketched.qs')

DefaultAssay(filtered) <- "sketch"

# perform clustering workflow
filtered <- FindVariableFeatures(filtered)
filtered <- ScaleData(filtered)
filtered <- RunPCA(filtered, assay = "sketch", reduction.name = "pca.sketch")
filtered <- FindNeighbors(filtered, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
filtered <- FindClusters(filtered, cluster.name = "seurat_cluster.sketched", resolution = .65)
filtered <- RunUMAP(filtered, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

Idents(filtered) <- "seurat_cluster.sketched"
DimPlot(filtered, reduction = "umap.sketch", label = F, cols = 'polychrome') + 
  ggtitle("Sketched clustering") + 
  theme(legend.position = "bottom")

DimPlot(filtered, reduction = "umap.sketch", label = F, cols = 'polychrome', split.by = 'sample', ncol = 4) + 
  ggtitle("Sketched clustering") + 
  theme(legend.position = "bottom")

meta <- filtered@meta.data
meta_sketch <- meta %>% filter(!is.na(seurat_cluster.sketched))

ggplot(meta_sketch, aes(y = general_code, fill = seurat_cluster.sketched)) + 
  geom_bar(position = 'fill') + scale_fill_cb_friendly()

filtered <- ProjectData(
  object = filtered,
  assay = "Spatial.016um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)


DefaultAssay(filtered) <- "Spatial.016um"
Idents(filtered) <- "seurat_cluster.projected"
qsave(filtered, 'data_processed/seurat_16um_sketched_projected.qs')

################# log normalized, all samples ##################################

# try with standard normalization
# result: clustering slightly less noisy

filtered <- NormalizeData(filtered)

filtered <- FindVariableFeatures(filtered)
filtered <- ScaleData(filtered)

filtered <- SketchData(
  object = filtered,
  assay = 'Spatial.016um',
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

qsave(filtered, 'data_processed/seurat_16um_norm_sketched.qs')

filtered <- qread('data_processed/seurat_16um_norm_sketched.qs')

DefaultAssay(filtered) <- "sketch"

# perform clustering workflow
filtered <- FindVariableFeatures(filtered)
filtered <- ScaleData(filtered)
filtered <- RunPCA(filtered, assay = "sketch", reduction.name = "pca.sketch")
filtered <- FindNeighbors(filtered, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
filtered <- FindClusters(filtered, cluster.name = "seurat_cluster.sketched", resolution = .65)
filtered <- RunUMAP(filtered, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

Idents(filtered) <- "seurat_cluster.sketched"
DimPlot(filtered, reduction = "umap.sketch", label = F, cols = 'polychrome') + 
  ggtitle("Sketched clustering") + 
  theme(legend.position = "bottom")

meta <- filtered@meta.data
meta_sketch <- meta %>% filter(!is.na(seurat_cluster.sketched))

ggplot(meta_sketch, aes(y = general_code, fill = seurat_cluster.sketched)) + 
  geom_bar(position = 'fill') + scale_fill_cb_friendly()

filtered <- ProjectData(
  object = filtered,
  assay = "Spatial.016um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)


DefaultAssay(filtered) <- "Spatial.016um"
Idents(filtered) <- "seurat_cluster.projected"
qsave(filtered, 'data_processed/seurat_16um_norm_sketched_projected.qs')

filtered <- qread('data_processed/seurat_16um_norm_sketched_projected.qs')

##############################log normalized, one sample #######################

# try to see if clustering is less noisy on only one sample as compared to all
# answer: not really

filtered_sample <- subset(filtered, subset = orig.ident == '2721_D7_Manual')

filtered_sample <- NormalizeData(filtered_sample)

filtered_sample <- FindVariableFeatures(filtered_sample)
filtered_sample <- ScaleData(filtered_sample)

filtered_sample <- SketchData(
  object = filtered_sample,
  assay = 'Spatial.016um',
  ncells = 10000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(filtered_sample) <- "sketch"

# perform clustering workflow
filtered_sample <- FindVariableFeatures(filtered_sample)
filtered_sample <- ScaleData(filtered_sample)
filtered_sample <- RunPCA(filtered_sample, assay = "sketch", reduction.name = "pca.sketch")
filtered_sample <- FindNeighbors(filtered_sample, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
filtered_sample <- FindClusters(filtered_sample, cluster.name = "seurat_cluster.sketched", resolution = .65)
filtered_sample <- RunUMAP(filtered_sample, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

Idents(filtered_sample) <- "seurat_cluster.sketched"
DimPlot(filtered_sample, reduction = "umap.sketch", label = F, cols = 'polychrome') + 
  ggtitle("Sketched clustering") + 
  theme(legend.position = "bottom")

filtered_sample <- ProjectData(
  object = filtered_sample,
  assay = "Spatial.016um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)


DefaultAssay(filtered_sample) <- "Spatial.016um"
Idents(filtered_sample) <- "seurat_cluster.projected"

color_pal = Seurat::DiscretePalette(n = length(unique(filtered_sample$seurat_cluster.projected)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(filtered_sample$seurat_cluster.projected))

SpatialDimPlot(filtered_sample, pt.size.factor = 5, cols = color_pal) #AGD13


################################################################################

# visualize results from various attempts above

color_pal = Seurat::DiscretePalette(n = length(unique(filtered$seurat_cluster.projected)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(filtered$seurat_cluster.projected))

SpatialDimPlot(filtered, images = 'X2721_D7_Manual.016um', pt.size.factor = 5, cols = color_pal) #AGD13
SpatialDimPlot(filtered, images = 'X3192Y_D2.016um', pt.size.factor = 5, cols = color_pal) #AGD3
SpatialDimPlot(filtered, images = 'Hct15HABA021_A5.016um', pt.size.factor = 5, cols = color_pal) #AGD8

SpatialDimPlot(filtered, images = 'X1013Y_A6.016um', pt.size.factor = 5, cols = color_pal) #AGD10
SpatialDimPlot(filtered, images = 'X845Y_A2.016um', pt.size.factor = 5, cols = color_pal) #AGD16



