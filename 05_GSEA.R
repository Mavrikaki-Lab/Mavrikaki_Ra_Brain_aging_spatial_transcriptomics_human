# 05_GSEA.R
# ------------------------------
# This script performs GSEA (Gene Set Enrichment Analysis) using fgsea
# for each spatial domain ("layer") from an Excel file with DEG results.
# Results are saved per layer and visualized across pathway collections.

# Load libraries
library(fgsea)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(tibble)
library(stringr)
library(viridis)
library(here)
library(readxl)

# Define paths
input_excel <- here("PathwayAnalysis", "AgeLayer_DEGs_Sex_Shrinkage.xlsx")  
output_base_dir <- here("PathwayAnalysis")

# Define GSEA collections to use
categories_df <- tibble::tibble(
  collection = c("C5", "C5", "C5", "C2", "C2", "H"),
  subcollection = c("GO:BP", "GO:CC", "GO:MF", "KEGG_LEGACY", "REACTOME", NA)
)

# GSEA function
run_gsea_analysis <- function(results, category = "C2", subcategory = NULL, species = "Homo sapiens",
                              minSize = 15, maxSize = 500, out_dir, scoreType = "std") {
  # Load gene sets
  pathways <- if (is.null(subcategory) || subcategory == "") {
    msigdbr(species = species, collection = category)
  } else {
    msigdbr(species = species, collection = category, subcollection = subcategory)
  }
  id <- paste0(category, if (!is.null(subcategory)) paste0("_", subcategory) else "")
  
  fgsea_sets <- pathways %>% split(x = .$gene_symbol, f = .$gs_name)
  
  # Prepare ranked list
  gsea_input <- results %>%
    arrange(desc(log2foldchange)) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    select(gene_name, log2foldchange)
  ranks <- deframe(gsea_input)
  
  # Run GSEA
  fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks, minSize = minSize,
                              maxSize = maxSize, eps = 0, scoreType = scoreType)
  
  # Process results
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidier <- fgseaResTidy %>%
    select(-leadingEdge, -ES, -log2err) %>%
    arrange(padj)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Save results
  write.csv(fgseaResTidier, file = file.path(out_dir, paste0("gsea_", id, ".csv")), row.names = FALSE)
  save(fgseaResTidy, file = file.path(out_dir, paste0("gsea_", id, ".RData")))
  
  # Save top/bottom pathways plot
  top_bottom <- bind_rows(
    fgseaResTidy %>% filter(padj < 0.05) %>% arrange(desc(NES)) %>% head(10),
    fgseaResTidy %>% filter(padj < 0.05) %>% arrange(NES) %>% head(10)
  )
  g <- ggplot(top_bottom, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = NES < 0)) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste("Top and Bottom GSEA Pathways -", id)) +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("Top_Bottom_GSEA_pathways_", id, ".pdf")),
         plot = g, width = 20, height = 12)
  
  return(fgseaResTidy)
}

# Run GSEA (layer/domain)
sheet_names <- readxl::excel_sheets(input_excel)
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

for (sheet_name in sheet_names) {
  message("Processing layer: ", sheet_name)
  
  ds <- read_excel(input_excel, sheet = sheet_name)
  colnames(ds) <- tolower(colnames(ds))
  
  if (!all(c("gene_name", "log2foldchange") %in% colnames(ds))) {
    message("Skipping ", sheet_name, ": missing required columns")
    next
  }
  
  ds <- ds %>%
    filter(!is.na(log2foldchange)) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  output_dir <- file.path(output_base_dir, sheet_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  gsea_results <- list()
  for (i in 1:nrow(categories_df)) {
    cat <- categories_df$collection[i]
    subcat <- categories_df$subcollection[i]
    
    message("  Running GSEA for: ", cat, " - ", subcat)
    gsea_results[[paste0(cat, "_", subcat)]] <- run_gsea_analysis(
      results = ds,
      category = cat,
      subcategory = ifelse(is.na(subcat), "", subcat),
      species = "Homo sapiens",
      minSize = 15,
      maxSize = 500,
      out_dir = output_dir
    )
  }
  
  saveRDS(gsea_results, file = file.path(output_dir, "gsea_all_results.rds"))
}

# Combined visualization per layer
layer_dirs <- list.dirs(output_base_dir, recursive = FALSE, full.names = TRUE)

for (layer_dir in layer_dirs) {
  rdata_files <- list.files(layer_dir, pattern = "\\.RData$", full.names = TRUE)
  if (length(rdata_files) == 0) next
  
  load_and_label <- function(file) {
    tryCatch({
      load(file)
      category_label <- gsub("gsea_|\\.RData", "", basename(file))
      fgseaResTidy %>%
        mutate(CAT = category_label) %>%
        arrange(desc(NES)) %>%
        slice_head(n = 10)
    }, error = function(e) {
      message("Error in ", file, ": ", e$message)
      return(NULL)
    })
  }
  
  top_pathways <- bind_rows(lapply(rdata_files, load_and_label)) %>%
    mutate(
      simplify_path = pathway %>%
        str_replace_all("GOBP_|GOCC_|GOMF_|HP_", "") %>%
        str_replace_all("_", " ") %>%
        str_replace_all("\\bREGULATION\\b", "REG") %>%
        str_replace_all("\\bINTRINSIC\\b", "INTR"),
      lp = -log10(pval)
    )
  
  p <- ggplot(top_pathways, aes(x = simplify_path, y = NES)) +
    geom_segment(aes(xend = simplify_path, yend = 0), color = "gray60", linewidth = 0.5) +
    geom_point(aes(size = size, color = lp), alpha = 0.9) +
    facet_wrap(vars(CAT), scales = "free_y", ncol = 2) +
    coord_flip() +
    scale_color_viridis_c(option = "magma", direction = -1) +
    scale_size_continuous(range = c(2, 6)) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.3, "in"),
      legend.position = "right"
    ) +
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      color = "-log10(p-value)",
      size = "Gene Set Size",
      title = paste("Top 10 GSEA Pathways -", basename(layer_dir))
    )
  
  ggsave(
    filename = file.path(layer_dir, paste0("Top10_GSEA_Pathways_", basename(layer_dir), ".pdf")),
    plot = p, width = 20, height = 12
  )
}
