#!/usr/bin/env Rscript

# Build Visium V2 pseudobulk counts by dominant cell2location type and BayesSpace domain.
# Counts are summed from whole spots; they are not purified cell-type RNA.
#
# Set data_dir below to the folder containing the two local input files.

# Load required libraries
library(SpatialExperiment)
library(Matrix)
library(here)

# Keep the input files outside GitHub and set data_dir to their local folder.
data_dir <- here("data", "cell2location_V2")
spe_file <- file.path(data_dir, "spe_filtered_final_BayesSpace8_Annotated.Rdata")
q05_file <- file.path(data_dir, "q05_21types_matched_to_annotated_SPE.csv")
output_dir <- here("results", "cell2location_V2_pseudobulk")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Read the annotated Visium V2 object.
loaded_objects <- new.env(parent = emptyenv())
load(spe_file, envir = loaded_objects)
spe <- loaded_objects$spe
stopifnot(inherits(spe, "SpatialExperiment"), "counts" %in% assayNames(spe))

metadata <- as.data.frame(colData(spe))
required_fields <- c("key", "sample_id", "group", "age", "sex", "BayesSpace", "cortical_layer")
stopifnot(all(required_fields %in% names(metadata)))
stopifnot(!anyNA(metadata[required_fields]), !anyDuplicated(metadata$key))

# Match q05 abundance rows to the spot order in the SpatialExperiment.
q05 <- read.csv(q05_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("key" %in% names(q05), setequal(as.character(metadata$key), as.character(q05$key)))
q05 <- q05[match(as.character(metadata$key), as.character(q05$key)), , drop = FALSE]
cell_types <- setdiff(names(q05), "key")
stopifnot(length(cell_types) == 21L, "Endothelial" %in% cell_types)

abundance <- as.matrix(q05[cell_types])
stopifnot(all(is.finite(abundance)), all(abundance >= 0), all(rowSums(abundance) > 0))

# Assign each spot to its highest-abundance cell type; ties use the first q05 column.
winner <- max.col(abundance, ties.method = "first")
tie_count <- sum(rowSums(abundance == abundance[cbind(seq_len(nrow(abundance)), winner)]) > 1L)
metadata$cell_type_c2l <- cell_types[winner]
metadata$pb_group <- paste(
  metadata$sample_id, metadata$cell_type_c2l, metadata$BayesSpace, sep = "_"
)

# Confirm the domain labels follow the Visium V2 layer mapping.
expected_layers <- c(
  Sp08D01 = "L1", Sp08D02 = "L2/3", Sp08D03 = "L3", Sp08D04 = "L3/4",
  Sp08D05 = "L5", Sp08D06 = "L6", Sp08D07 = "WM1", Sp08D08 = "WM2"
)
domain <- as.character(metadata$BayesSpace)
stopifnot(all(domain %in% names(expected_layers)))
stopifnot(all(as.character(metadata$cortical_layer) == unname(expected_layers[domain])))

# Sample-level annotations should be constant across spots from the same sample.
for (field in intersect(
  c("group", "age", "sex", "subject", "tissue_bank", "processing_batch"),
  names(metadata)
)) {
  values_by_sample <- split(as.character(metadata[[field]]), metadata$sample_id)
  stopifnot(!anyNA(metadata[[field]]))
  stopifnot(all(vapply(values_by_sample, function(x) length(unique(x)) == 1L, logical(1))))
}

# Sum spot counts for each sample, assigned cell type, and domain.
counts <- assay(spe, "counts")
genes <- rownames(spe)
stopifnot(!is.null(genes), !anyNA(genes), !anyDuplicated(genes))

# Check counts in blocks to avoid making a second dense copy of the full assay.
gene_blocks <- split(seq_len(nrow(counts)), ceiling(seq_len(nrow(counts)) / 500L))
for (genes_in_block in gene_blocks) {
  count_block <- as.matrix(counts[genes_in_block, , drop = FALSE])
  stopifnot(all(is.finite(count_block)), all(count_block >= 0),
            all(count_block == round(count_block)))
}

pseudobulk_groups <- sort(unique(metadata$pb_group))
spot_to_group <- sparseMatrix(
  i = seq_len(nrow(metadata)),
  j = match(metadata$pb_group, pseudobulk_groups),
  x = 1,
  dims = c(nrow(metadata), length(pseudobulk_groups))
)
message("Aggregating ", ncol(counts), " spots into ", length(pseudobulk_groups), " pseudobulks...")
counts_pb <- matrix(
  0,
  nrow = nrow(counts),
  ncol = length(pseudobulk_groups),
  dimnames = list(genes, pseudobulk_groups)
)
for (genes_in_block in gene_blocks) {
  counts_pb[genes_in_block, ] <- as.matrix(counts[genes_in_block, , drop = FALSE] %*% spot_to_group)
}
stopifnot(all(rowSums(counts_pb) == rowSums(counts)), all(colSums(counts_pb) > 0))

metadata_columns <- intersect(
  c("pb_group", "sample_id", "subject", "group", "age", "sex", "tissue_bank",
    "processing_batch", "cell_type_c2l", "BayesSpace", "cortical_layer",
    "cortical_layer_merged", "maynard_naming"),
  names(metadata)
)
meta_pb2 <- metadata[match(pseudobulk_groups, metadata$pb_group), metadata_columns, drop = FALSE]
rownames(meta_pb2) <- pseudobulk_groups
meta_pb2$age_layer <- paste(meta_pb2$group, meta_pb2$BayesSpace, sep = "_")
meta_pb2$cell_state <- paste(meta_pb2$cell_type_c2l, meta_pb2$BayesSpace, sep = "_")
meta_pb2$n_spots <- tabulate(match(metadata$pb_group, pseudobulk_groups), nbins = length(pseudobulk_groups))
stopifnot(identical(colnames(counts_pb), rownames(meta_pb2)))
stopifnot(sum(meta_pb2$n_spots) == ncol(spe))

saveRDS(counts_pb, file.path(output_dir, "counts_pb.rds"))
saveRDS(meta_pb2, file.path(output_dir, "meta_pb2.rds"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Pseudobulk files written to: ", output_dir)






