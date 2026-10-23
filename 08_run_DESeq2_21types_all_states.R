#!/usr/bin/env Rscript
# Reproduce the Visium V2 cell-type/domain DESeq2 analysis.
# Whole-spot RNA was grouped by dominant cell2location type and BayesSpace domain.

# Load required libraries
library(here)
library(DESeq2)
library(ashr)

# Set input and output directories
input_dir <- here("cell2location_DE_reproduction", "data", "pseudobulk_21types")
output_dir <- here("cell2location_DE_reproduction", "results", "DESeq2_21types_all_states")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read pseudobulk counts and metadata
counts_pb <- readRDS(file.path(input_dir, "counts_pb.rds"))
meta_pb <- as.data.frame(readRDS(file.path(input_dir, "meta_pb2.rds")))

meta_pb$cell_state <- paste(meta_pb$cell_type_c2l, meta_pb$BayesSpace, sep = "_")
rownames(meta_pb) <- meta_pb$pb_group

meta_pb <- meta_pb[match(colnames(counts_pb), meta_pb$pb_group), , drop = FALSE]
# Make a comparison plan for all cell types and domains.
domains <- c(
  Sp08D01 = "L1", Sp08D02 = "L2/3", Sp08D03 = "L3", Sp08D04 = "L3/4",
  Sp08D05 = "L5", Sp08D06 = "L6", Sp08D07 = "WM1", Sp08D08 = "WM2"
)
cell_types <- sort(unique(meta_pb$cell_type_c2l))
plan <- expand.grid(
  cell_type = cell_types,
  domain = names(domains),
  stringsAsFactors = FALSE
)
plan$cell_state <- paste(plan$cell_type, plan$domain, sep = "_")
plan$layer <- unname(domains[plan$domain])
plan$nYoung <- 0L
plan$nOld <- 0L
plan$nMiddle <- 0L
plan$min_spots <- NA_real_
plan$pseudobulks_under10spots <- 0L
  c("Astro_FB", "Astro_PP", "Exc_L2-3", "Endothelial"), "Sp08D02", sep = "_"
)
plan$status <- "pending"
plan$reason <- ""
plan$genes <- NA_integer_
plan$genes_with_padj <- NA_integer_
plan$n_up <- NA_integer_
plan$n_down <- NA_integer_

for (i in seq_len(nrow(plan))) {
  state_metadata <- meta_pb[meta_pb$cell_state == plan$cell_state[i], , drop = FALSE]
  age_metadata <- state_metadata[state_metadata$group %in% c("Young", "Old"), , drop = FALSE]
  n_young <- sum(age_metadata$group == "Young")
  n_old <- sum(age_metadata$group == "Old")

  plan$nYoung[i] <- n_young
  plan$nOld[i] <- n_old
  plan$nMiddle[i] <- sum(state_metadata$group == "Middle")
  plan$min_spots[i] <- if (nrow(age_metadata) > 0L) min(age_metadata$n_spots) else NA_real_
  plan$pseudobulks_under10spots[i] <- sum(age_metadata$n_spots < 10)

  if (n_young < 2L || n_old < 2L) {
    plan$status[i] <- "skipped"
    plan$reason[i] <- "Fewer than two Young or Old pseudobulks"
  }
}

write.csv(plan, file.path(output_dir, "state_status.csv"), row.names = FALSE)

# Use the same global count filter as the validated analysis.
keep_samples <- meta_pb$group %in% c("Young", "Old")
counts_pb <- counts_pb[, keep_samples, drop = FALSE]
meta_pb <- meta_pb[keep_samples, , drop = FALSE]
gene_keep <- rowSums(counts_pb >= 10) >= 3
counts_pb <- counts_pb[gene_keep, , drop = FALSE]
results_by_state <- list()

for (i in which(plan$status == "pending")) {
  state <- plan$cell_state[i]
  message("Running ", state, " (Young n = ", plan$nYoung[i], ", Old n = ", plan$nOld[i], ")")

  tryCatch({
    state_samples <- meta_pb$cell_state == state
    state_metadata <- meta_pb[state_samples, , drop = FALSE]
    state_counts <- counts_pb[, state_samples, drop = FALSE]
    state_metadata$age_group <- factor(state_metadata$group, levels = c("Young", "Old"))

    state_counts <- state_counts[rowSums(state_counts) > 0, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(
      countData = round(state_counts),
      colData = state_metadata,
      design = ~ age_group
    )
    dds <- DESeq(dds, parallel = FALSE)
    shrunk_results <- lfcShrink(
      dds,
      contrast = c("age_group", "Old", "Young"),
      type = "ashr"
    )

    result_table <- as.data.frame(shrunk_results)
    result_table$gene_id <- rownames(result_table)
    result_table$cell_state <- state
    result_table$cell_type <- plan$cell_type[i]
    result_table$BayesSpace <- plan$domain[i]
    result_table$cortical_layer <- plan$layer[i]
    result_table$contrast <- "Old vs Young"
    results_by_state[[state]] <- result_table

    plan$genes[i] <- nrow(result_table)
    plan$genes_with_padj[i] <- sum(!is.na(result_table$padj))
    plan$n_up[i] <- sum(result_table$padj < 0.05 & result_table$log2FoldChange > 0.25, na.rm = TRUE)
    plan$n_down[i] <- sum(result_table$padj < 0.05 & result_table$log2FoldChange < -0.25, na.rm = TRUE)
    plan$status[i] <- "completed"

    rm(dds, raw_results, shrunk_results, state_counts)
    gc(verbose = FALSE)
  }, error = function(e) {
    plan$status[i] <<- "failed"
    plan$reason[i] <<- conditionMessage(e)
    message("Failed ", state, ": ", conditionMessage(e))
  })

  write.csv(plan, file.path(output_dir, "state_status.csv"), row.names = FALSE)
}

if (length(results_by_state) > 0L) {
  all_results <- do.call(rbind, results_by_state)
  rownames(all_results) <- NULL
  saveRDS(all_results, file.path(output_dir, "all_states_ashr_results.rds"))
}

write.csv(plan, file.path(output_dir, "state_status.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message(
  "Completed: ", sum(plan$status == "completed"),
  "; skipped: ", sum(plan$status == "skipped"),
  "; failed: ", sum(plan$status == "failed"),
  ". Results: ", output_dir
)
if (any(plan$status == "failed")) stop("Some cell states failed; see state_status.csv.")
message("All-state ashr analysis complete.")
