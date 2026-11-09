library(tidyverse)
library(Seurat)
library(qs)
library(BPCells)

# this script needs ~150GB of memory to run


####################### normalized, sketched, projected ########################

norm_sketched_proj <- qread('data_processed/seurat_16um_norm_sketched_projected.qs')

write_matrix_dir(mat = norm_sketched_proj[["Spatial.016um"]]$counts, dir = 'data_processed/seurat_16um_norm_sketched_projected_matrix_counts/')
write_matrix_dir(mat = norm_sketched_proj[["Spatial.016um"]]$data, dir = 'data_processed/seurat_16um_norm_sketched_projected_matrix_data/')

counts.mat <- open_matrix_dir(dir = 'data_processed/seurat_16um_norm_sketched_projected_matrix_counts/')
data.mat <- open_matrix_dir(dir = 'data_processed/seurat_16um_norm_sketched_projected_matrix_data/')

norm_sketched_proj[["Spatial.016um"]]$counts <- counts.mat
norm_sketched_proj[["Spatial.016um"]]$data <- data.mat
norm_sketched_proj[["Spatial.016um"]]$scale.data <- NULL

qsave(norm_sketched_proj, 'data_processed/seurat_16um_norm_sketched_projected_ondisk.qs')

rm(norm_sketched_proj)
gc()

############### non-normalized, sketched, projected #############################

sketched_proj <- qread('data_processed/seurat_16um_sketched_projected.qs')

write_matrix_dir(mat = sketched_proj[["Spatial.016um"]]$counts, dir = 'data_processed/seurat_16um_sketched_projected_matrix_counts/')
write_matrix_dir(mat = sketched_proj[["Spatial.016um"]]$data, dir = 'data_processed/seurat_16um_sketched_projected_matrix_data/')

counts.mat <- open_matrix_dir(dir = 'data_processed/seurat_16um_sketched_projected_matrix_counts/')
data.mat <- open_matrix_dir(dir = 'data_processed/seurat_16um_sketched_projected_matrix_data/')

sketched_proj[["Spatial.016um"]]$counts <- counts.mat
sketched_proj[["Spatial.016um"]]$data <- data.mat
sketched_proj[["Spatial.016um"]]$scale.data <- NULL

qsave(sketched_proj, 'data_processed/seurat_16um_sketched_projected_ondisk.qs')
