# library(sceasy)
# library(reticulate)
# library(Seurat)
# 
# # to make conda environment, used terminal interactive session in FAS
# 
# use_condaenv('~/.conda/envs/sceasy')
# 
# # doesn't work,
# # Error in py_module_import(module, convert = convert) : 
# # ImportError: /usr/lib/x86_64-linux-gnu/libssl.so.3: version `OPENSSL_3.2.0' not found (required by /n/home12/abartlett04/.conda/envs/sceasy/lib/python3.9/site-packages/h5py/../../.././libcurl.so.4)
# 
# sceasy::convertFormat('data/sce_binder.h5ad.gz', from="anndata", to="seurat",
#                       outFile='sce_binder.rds')


################################################################################

# fails in Convert command, Error: Not a sparse matrix

# library(SeuratDisk)
# 
# Convert('data/sce_binder.h5ad', dest = "h5seurat", overwrite = TRUE)
# sce_binder <- LoadH5Seurat('data/sce_binder.h5seurat')

################################################################################

library(Matrix)
library(qs)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(Seurat)

################ using logcounts of ref ########################################

# uses output of prep_rctd_ref.ipynb to complete conversion from anndata to seurat
save_dir <- 'data/sce_binder_for_R/'

counts <- readMM(paste0(save_dir,'/counts.mtx'))
dim(counts)
cellMeta <- read.csv(paste0(save_dir,'/counts_cellMeta.csv'))
head(cellMeta)
geneMeta <- read.csv(paste0(save_dir,'/counts_geneMeta.csv'))
dim(geneMeta)
head(geneMeta)

# Set the rownames and colnames
rownames(counts) <- cellMeta$Barcode
colnames(counts) <- geneMeta$GeneName

# filter counts to have only genes with names (to match our dataset)
geneMeta_new <- AnnotationDbi::select(org.Hs.eg.db, geneMeta$GeneName, "ENSEMBL", columns = "SYMBOL") %>%
  filter(!duplicated(ENSEMBL)) %>% 
  mutate(keep = !is.na(SYMBOL) & !duplicated(SYMBOL)) 

counts_new <- counts[,geneMeta_new$keep]
colnames(counts_new) <- geneMeta_new$SYMBO[geneMeta_new$keep]

# create seurat object, set metadata
seo <- CreateSeuratObject(counts = t(counts_new), project = "min", min.cells = 3, min.features = 200)
seo@meta.data <- cbind(cellMeta,seo@meta.data)
rownames(seo@meta.data) <- colnames(seo)

qsave(seo, 'data/sce_binder.qs')

################ using raw counts of ref ########################################

# uses output of prep_rctd_ref.ipynb to complete conversion from anndata to seurat
save_dir <- 'data/sce_binder_for_R_raw/'

counts <- readMM(paste0(save_dir,'/counts.mtx'))
dim(counts)
cellMeta <- read.csv(paste0(save_dir,'/counts_cellMeta.csv'))
head(cellMeta)
geneMeta <- read.csv(paste0(save_dir,'/counts_geneMeta.csv'))
dim(geneMeta)
head(geneMeta)

# Set the rownames and colnames
rownames(counts) <- cellMeta$Barcode
colnames(counts) <- geneMeta$GeneName

# filter counts to have only genes with names (to match our dataset)
geneMeta_new <- AnnotationDbi::select(org.Hs.eg.db, geneMeta$GeneName, "ENSEMBL", columns = "SYMBOL") %>%
  filter(!duplicated(ENSEMBL)) %>% 
  mutate(keep = !is.na(SYMBOL) & !duplicated(SYMBOL)) 

counts_new <- counts[,geneMeta_new$keep]
colnames(counts_new) <- geneMeta_new$SYMBO[geneMeta_new$keep]

# create seurat object, set metadata
seo <- CreateSeuratObject(counts = t(counts_new), project = "min", min.cells = 3, min.features = 200)
seo@meta.data <- cbind(cellMeta,seo@meta.data)
rownames(seo@meta.data) <- colnames(seo)

qsave(seo, 'data/sce_binder_raw.qs')
