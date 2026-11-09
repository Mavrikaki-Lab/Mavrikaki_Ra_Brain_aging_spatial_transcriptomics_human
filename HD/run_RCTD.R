library(spacexr)
library(Seurat)
library(BPCells)
library(tidyverse)
library(qs)
options(future.globals.maxSize= 20000000000)

##################### run RCTD #################################################

# ref_fn <- 'data/sce_binder.qs' # using logcounts of ref
ref_fn <- 'data/sce_binder_raw.qs' # using raw counts of ref

seurat_norm_clust <- qread('data_processed/seurat_16um_norm_sketched_projected_ondisk.qs')
seurat_ref <- qread(ref_fn) 

# create the RCTD query object

coords <- lapply(names(seurat_norm_clust@images), function(image_name) {
  GetTissueCoordinates(seurat_norm_clust, image = image_name)
}) %>% bind_rows()

seurat_norm_clust <- subset(seurat_norm_clust, cells = rownames(coords))
counts <- seurat_norm_clust[["sketch"]]$counts
cells <- colnames(counts)

coords <- coords[cells, 1:2]
query <- SpatialRNA(coords, counts, colSums(counts))

# create the RCTD reference object
Idents(seurat_ref) <- "celltypes_final"
counts_ref <- seurat_ref[["RNA"]]$counts
cluster <- as.factor(seurat_ref$celltypes_final)
nUMI <- seurat_ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

reference <- Reference(counts_ref, cluster, nUMI)

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 64)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

if(grepl('raw', ref_fn)){
  qsave(RCTD, 'data_processed/RCTD_raw.qs')
} else {
  qsave(RCTD, 'data_processed/RCTD.qs')
}

RCTD_meta <- RCTD@results$results_df
RCTD_meta$barcode <- rownames(RCTD_meta)



seurat_norm_clust <- AddMetaData(seurat_norm_clust, metadata = RCTD_meta)
seurat_norm_clust$first_type <- as.character(seurat_norm_clust$first_type)
seurat_norm_clust$first_type[is.na(seurat_norm_clust$first_type)] <- "Unknown"


seurat_norm_clust <- ProjectData(
  object = seurat_norm_clust,
  assay = "Spatial.016um",
  full.reduction = "full.pca.sketch.rctd",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type"),
  recompute.weights = T
)

color_pal = Seurat::DiscretePalette(n = length(unique(seurat_norm_clust$full_first_type)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(seurat_norm_clust$full_first_type))
SpatialDimPlot(seurat_norm_clust, group.by = 'full_first_type', 
               pt.size.factor = 7, cols = color_pal, ncol = 5) 
if(grepl('raw', ref_fn)){
  qsave(seurat_norm_clust, 'data_processed/seurat_16um_norm_sketched_projected_ondisk_rctd_raw.qs')
} else{
qsave(seurat_norm_clust, 'data_processed/seurat_16um_norm_sketched_projected_ondisk_rctd.qs')
}

########################### compare RCTD results ###############################

# takeaway: raw seems better, not surprising
raw <- qread('data_processed/seurat_16um_norm_sketched_projected_ondisk_rctd_raw.qs')
log <- qread('data_processed/seurat_16um_norm_sketched_projected_ondisk_rctd.qs')

coords <- lapply(names(raw@images), function(image_name) {
  GetTissueCoordinates(raw, image = image_name)
}) %>% bind_rows()

raw <- AddMetaData(raw, coords)

meta <- raw@meta.data %>% as.data.frame %>% 
  left_join(log@meta.data %>% as.data.frame %>% dplyr::select(barcode, full_first_type_log = full_first_type)) %>%
  pivot_longer(c(full_first_type, full_first_type_log), names_to = 'normalization', values_to = 'cell_type')

color_pal = Seurat::DiscretePalette(n = length(unique(meta$cell_type)),
                                    palette = "polychrome")
names(color_pal) <- sort(unique(meta$cell_type))

for (sample_name in sort(unique(meta$general_code))){
  p <- ggplot(meta %>% filter(general_code == sample_name), aes(x = x, y = y, color = cell_type)) + 
    facet_wrap(~normalization) +
    geom_point(size = .15, shape = 'square') + 
    scale_color_manual(values = color_pal) + 
    guides(color = guide_legend(override.aes = list(size=3, alpha=1), ncol = 2)) + 
    ggtitle(sample_name)
  print(p)
}

meta_raw <- meta %>% filter(normalization == 'full_first_type')

for (sample_name in sort(unique(meta_raw$general_code))){
  p <- ggplot(meta_raw %>% filter(general_code == sample_name), aes(x = x, y = y, color = cell_type)) + 
    facet_wrap(~cell_type, ncol = 5) +
    geom_point(size = .15, shape = 'square') + 
    scale_color_manual(values = color_pal) + 
    guides(color = guide_legend(override.aes = list(size=3, alpha=1), ncol = 2)) + 
    ggtitle(sample_name)
  print(p)
}
