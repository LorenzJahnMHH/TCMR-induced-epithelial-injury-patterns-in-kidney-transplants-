# Suppl. Fig. 3
Idents(Ktx_data_mouse) <- "ID_to_plot"

ordered_ids <- c(
  "BALB/c to BALB/c 1", "BALB/c to BALB/c 2",  "C57BL/6 to C57BL/6 1", "BALB/c to C57BL/6 1",
  "BALB/c to C57BL/6 2", "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3"
)

Ktx_data_mouse$ID_to_plot <- factor(Ktx_data_mouse$ID_to_plot, levels = ordered_ids)

sample_colors <- c(
  "BALB/c to BALB/c 1" = "grey85",  "BALB/c to BALB/c 2" = "grey85",  "C57BL/6 to C57BL/6 1" = "grey85",
  "BALB/c to C57BL/6 1" = "#B4EEB4", "BALB/c to C57BL/6 2" = "#B4EEB4",
  "C57BL/6 to BALB/c 1" = "#BCD2EE", "C57BL/6 to BALB/c 2" = "#BCD2EE", "C57BL/6 to BALB/c 3" = "#BCD2EE"
)

features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")

for (feature in features) {
  p <- VlnPlot(
    Ktx_data_mouse, 
    features = c(feature), 
    pt.size = 0
  ) + 
    scale_fill_manual(values = sample_colors) +
    theme(
      axis.text.x = element_text(size = rel(0.8))
    ) +
    NoLegend()
  
  print(p)
}

library(ggplot2)

cell_counts <- Ktx_data_mouse@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(cell_count = n())

ggplot(cell_counts, aes(x = ID_to_plot, y = cell_count)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    axis.text.y = element_text(size = 10)
  ) +
  ylab("")

# Suppl. Fig 6A
balbc_dir <- "D:/R_Ordner/KTx Daten/Pathway_enrichment/Balbc_down_pct_01/"
blck6_dir <- "D:/R_Ordner/KTx Daten/Pathway_enrichment/Blck6_down_pct_01/"

process_enrichr_file <- function(file_path) {
  data <- read.delim(file_path, header = TRUE, sep = "\t", fill = TRUE, quote = "", check.names = FALSE)
  clean_data <- data %>%
    dplyr::select(Term, `Adjusted P-value`) %>%
    dplyr::rename(Pathway = Term, AdjustedPValue = `Adjusted P-value`) %>%
    dplyr::mutate(
      Pathway = str_trim(Pathway),
      Pathway = str_replace_all(Pathway, "\\s+", " "),
      Pathway = sub("^Pp", "P", Pathway),
      AdjustedPValue = as.numeric(AdjustedPValue)
    )
  return(clean_data)
}

celltypes <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "IntC", "Uro", "PEC")

all_pathway_data <- tibble(
  Pathway = character(),
  AdjustedPValue = numeric(),
  CellType = character(),
  Strain = character()
)

for (ct in celltypes) {
  balbc_file <- paste0(balbc_dir, "Balbc_", ct, "_pathways_raw.txt")
  blck6_file <- paste0(blck6_dir, "Blck6_", ct, "_pathways_raw.txt")
  
  if (file.exists(balbc_file)) {
    balbc_data <- process_enrichr_file(balbc_file) %>%
      mutate(CellType = ct, Strain = "Balbc")
    all_pathway_data <- bind_rows(all_pathway_data, balbc_data)
  }
  
  if (file.exists(blck6_file)) {
    blck6_data <- process_enrichr_file(blck6_file) %>%
      mutate(CellType = ct, Strain = "Blck6")
    all_pathway_data <- bind_rows(all_pathway_data, blck6_data)
  }
}

filtered_data <- all_pathway_data %>%
  filter(!is.na(AdjustedPValue), AdjustedPValue < 0.05) %>%
  mutate(
    log10_fdr = -log10(AdjustedPValue),
    CellType = factor(CellType, levels = celltypes),
    Strain   = factor(Strain, levels = c("Balbc", "Blck6"))
  )

reshaped_data <- filtered_data %>%
  unite("CellType_Strain", CellType, Strain, sep = "_", remove = FALSE) %>%
  pivot_wider(
    names_from  = CellType_Strain,
    values_from = log10_fdr,
    values_fill = 0
  )

pathway_counts <- filtered_data %>%
  group_by(Pathway) %>%
  summarize(celltype_count = sum(log10_fdr > 0, na.rm = TRUE), .groups = "drop") %>%
  arrange(celltype_count)

top_pathways <- filtered_data %>%
  group_by(Pathway) %>%
  summarize(min_adj_pval = min(AdjustedPValue, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_adj_pval) %>%
  slice_head(n = 25)

filtered_top_data <- filtered_data %>%
  filter(Pathway %in% top_pathways$Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = pathway_counts$Pathway))

dotplot_sorted_reversed <- ggplot(filtered_top_data, aes(x = Strain,
                                                         y = Pathway,
                                                         size = pmin(log10_fdr, 3),
                                                         fill = Strain)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) +
  scale_size_continuous(range = c(0, 6), breaks = c(0, 1, 2, 3), limits = c(0, 3)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
    strip.text = element_blank(),
    panel.spacing.x = unit(2, "lines"),
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  ) +
  facet_grid(~CellType, scales = "free_x", space = "free_x")

print(dotplot_sorted_reversed)

dotplot_sorted <- ggplot(filtered_top_data, aes(x = Strain,
                                                y = Pathway,
                                                size = pmin(log10_fdr, 3),
                                                fill = Strain)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) +
  scale_size_continuous(range = c(0, 6), breaks = c(0, 1, 2, 3), limits = c(0, 3)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing.x = unit(2, "lines"),
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  ) +
  labs(
    x = "Strain",
    y = "Pathway",
    size = "-log10(FDR)",
    fill = "Strain"
  ) +
  facet_grid(~CellType, scales = "free_x", space = "free_x")

print(dotplot_sorted)

# Suppl. Fig. 6 B
setwd(".../Blck6_Balbc_DGE_lists")


library(Seurat)
library(tidyverse)
library(pheatmap)

Balbc_up_Allo_vs_Syn <- readRDS("Balbc_up_Allo_vs_Syn.rds")
Balbc_down_Allo_vs_Syn <- readRDS("Balbc_down_Allo_vs_Syn.rds")
Blck6_up_Allo_vs_Syn <- readRDS("Blck6_up_Allo_vs_Syn.rds")
Blck6_down_Allo_vs_Syn <- readRDS("Blck6_down_Allo_vs_Syn.rds")

remove_celltype <- "Leuko"
Balbc_up_Allo_vs_Syn <- Balbc_up_Allo_vs_Syn[!names(Balbc_up_Allo_vs_Syn) %in% remove_celltype]
Balbc_down_Allo_vs_Syn <- Balbc_down_Allo_vs_Syn[!names(Balbc_down_Allo_vs_Syn) %in% remove_celltype]
Blck6_up_Allo_vs_Syn <- Blck6_up_Allo_vs_Syn[!names(Blck6_up_Allo_vs_Syn) %in% remove_celltype]
Blck6_down_Allo_vs_Syn <- Blck6_down_Allo_vs_Syn[!names(Blck6_down_Allo_vs_Syn) %in% remove_celltype]

unique_genes_up <- sort(unique(c(
  unlist(lapply(Blck6_up_Allo_vs_Syn, rownames)), 
  unlist(lapply(Balbc_up_Allo_vs_Syn, rownames))
)))

unique_genes_down <- sort(unique(c(
  unlist(lapply(Blck6_down_Allo_vs_Syn, rownames)), 
  unlist(lapply(Balbc_down_Allo_vs_Syn, rownames))
)))

integrated_matrix_up <- data.frame(gene = unique_genes_up)
integrated_matrix_down <- data.frame(gene = unique_genes_down)

for (cell_type in names(Blck6_up_Allo_vs_Syn)) {
  integrated_matrix_up[paste0("Up_", cell_type, "_Blck6")] <- NA
  integrated_matrix_down[paste0("Down_", cell_type, "_Blck6")] <- NA
  integrated_matrix_up[paste0("Up_", cell_type, "_Balbc")] <- NA
  integrated_matrix_down[paste0("Down_", cell_type, "_Balbc")] <- NA
}

for (cell_type in names(Blck6_up_Allo_vs_Syn)) {
  if (nrow(Blck6_up_Allo_vs_Syn[[cell_type]]) > 0) {
    match_idx <- match(rownames(Blck6_up_Allo_vs_Syn[[cell_type]]), integrated_matrix_up$gene)
    valid_idx <- !is.na(match_idx)
    integrated_matrix_up[match_idx[valid_idx], paste0("Up_", cell_type, "_Blck6")] <- Blck6_up_Allo_vs_Syn[[cell_type]]$avg_log2FC[valid_idx]
  }
  if (nrow(Blck6_down_Allo_vs_Syn[[cell_type]]) > 0) {
    match_idx <- match(rownames(Blck6_down_Allo_vs_Syn[[cell_type]]), integrated_matrix_down$gene)
    valid_idx <- !is.na(match_idx)
    integrated_matrix_down[match_idx[valid_idx], paste0("Down_", cell_type, "_Blck6")] <- Blck6_down_Allo_vs_Syn[[cell_type]]$avg_log2FC[valid_idx]
  }
}

for (cell_type in names(Balbc_up_Allo_vs_Syn)) {
  if (nrow(Balbc_up_Allo_vs_Syn[[cell_type]]) > 0) {
    match_idx <- match(rownames(Balbc_up_Allo_vs_Syn[[cell_type]]), integrated_matrix_up$gene)
    valid_idx <- !is.na(match_idx)
    integrated_matrix_up[match_idx[valid_idx], paste0("Up_", cell_type, "_Balbc")] <- Balbc_up_Allo_vs_Syn[[cell_type]]$avg_log2FC[valid_idx]
  }
  if (nrow(Balbc_down_Allo_vs_Syn[[cell_type]]) > 0) {
    match_idx <- match(rownames(Balbc_down_Allo_vs_Syn[[cell_type]]), integrated_matrix_down$gene)
    valid_idx <- !is.na(match_idx)
    integrated_matrix_down[match_idx[valid_idx], paste0("Down_", cell_type, "_Balbc")] <- Balbc_down_Allo_vs_Syn[[cell_type]]$avg_log2FC[valid_idx]
  }
}

integrated_matrix_up[is.na(integrated_matrix_up)] <- 0
integrated_matrix_down[is.na(integrated_matrix_down)] <- 0

rownames(integrated_matrix_up) <- integrated_matrix_up$gene
integrated_matrix_up$gene <- NULL

rownames(integrated_matrix_down) <- integrated_matrix_down$gene
integrated_matrix_down$gene <- NULL

quantile_up <- quantile(as.numeric(integrated_matrix_up[integrated_matrix_up > 0]), probs = 0.7, na.rm = TRUE)
quantile_down <- quantile(as.numeric(integrated_matrix_down[integrated_matrix_down < 0]), probs = 0.3, na.rm = TRUE)



breaks_up <- seq(0, quantile_up, length.out = 256)
color_palette_up <- colorRampPalette(c("beige", "darkorange"))(255)

breaks_down <- seq(quantile_down, 0, length.out = 256)
color_palette_down <- colorRampPalette(c("royalblue", "beige"))(255)


column_order_up <- c()
column_order_down <- c()

for (cell_type in names(Blck6_up_Allo_vs_Syn)) {
  column_order_up <- c(column_order_up, paste0("Up_", cell_type, "_Balbc"), paste0("Up_", cell_type, "_Blck6"))
}

for (cell_type in names(Blck6_down_Allo_vs_Syn)) {
  column_order_down <- c(column_order_down, paste0("Down_", cell_type, "_Balbc"), paste0("Down_", cell_type, "_Blck6"))
}

cell_types <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "IntC", "Uro", "PEC")
cell_type_labels_up <- cell_types 
cell_type_labels_down <- cell_types  

pheatmap(
  integrated_matrix_up[, column_order_up, drop = FALSE],
  scale = "none",
  clustering_distance_rows = "correlation", 
  cluster_cols = FALSE,
  color = color_palette_up,
  na_col = "white",
  breaks = breaks_up,
  show_rownames = TRUE, 
  show_colnames = TRUE,
  labels_col = cell_type_labels_up
)

pheatmap(
  integrated_matrix_down[, column_order_down, drop = FALSE],
  scale = "none",
  clustering_distance_rows = "correlation", 
  cluster_cols = FALSE,
  color = color_palette_down,
  na_col = "white",
  breaks = breaks_down,
  show_rownames = TRUE, 
  show_colnames = TRUE,
  labels_col = cell_type_labels_down
)

pheatmap(
  integrated_matrix_up[, column_order_up, drop = FALSE],
  scale = "none",
  clustering_distance_rows = "correlation", 
  cluster_cols = FALSE,
  color = color_palette_up,
  na_col = "white",
  breaks = breaks_up,
  show_rownames = FALSE, 
  show_colnames = FALSE,
  labels_col = NULL,
  legend = FALSE
)

pheatmap(
  integrated_matrix_down[, column_order_down, drop = FALSE],
  scale = "none",
  clustering_distance_rows = "correlation", 
  cluster_cols = FALSE,
  color = color_palette_down,
  na_col = "white",
  breaks = breaks_down,
  show_rownames = FALSE, 
  show_colnames = FALSE,
  labels_col = NULL,
  legend = FALSE
)


# Suppl. Fig. 7 A+B
library(WGCNA)
library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(forcats)
library(tidyr)
library(stringr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

Blck6_up_Allo_vs_Syn <- readRDS("Blck6_up_Allo_vs_Syn.rds")
Balbc_up_Allo_vs_Syn <- readRDS("Balbc_up_Allo_vs_Syn.rds")

ct_keep <- c("Podo", "PT", "tL", "TAL", "CNT", "DCT",
             "CD-IC-A", "CD-IC-B", "CD-PC", "PEC")

lfc_matrix_from_deglist <- function(deg_list, keep_ct) {
  common_genes <- Reduce(intersect, lapply(deg_list, rownames))
  deg_list <- lapply(deg_list, function(df) df[common_genes, , drop = FALSE])
  mat <- do.call(cbind, lapply(deg_list[keep_ct], function(df) df$avg_log2FC))
  colnames(mat) <- keep_ct
  rownames(mat) <- common_genes
  return(mat)
}

pretty_pathway_name <- function(description) {
  if (startsWith(description, "REACTOME_")) {
    rest <- gsub("REACTOME_", "", description)
    return(paste(stringr::str_to_title(gsub("_", " ", rest)), "- RC"))
  } else if (startsWith(description, "WP_")) {
    rest <- gsub("WP_", "", description)
    return(paste(stringr::str_to_title(gsub("_", " ", rest)), "- WP"))
  } else if (startsWith(description, "BIOCARTA_")) {
    rest <- gsub("BIOCARTA_", "", description)
    return(paste(stringr::str_to_title(gsub("_", " ", rest)), "- BC"))
  } else if (startsWith(description, "HALLMARK_")) {
    rest <- gsub("HALLMARK_", "", description)
    return(paste(stringr::str_to_title(gsub("_", " ", rest)), "- HM"))
  } else {
    return(description)
  }
}

# adjust for pathway-databases and q-values for 6C, see methods
run_wgcna_from_lfc <- function(lfc_mat, prefix,
                               minModuleSize = 20, deepSplit = 3, mergeCutHeight = 0.20) {
  datExpr <- t(lfc_mat)
  datExpr <- datExpr[, apply(datExpr, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  
  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
  
  pdf(paste0(prefix, "_pickSoftThreshold.pdf"), 7, 4)
  par(mfrow = c(1, 2)); cex1 <- 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels = powers, cex = cex1, col = "red")
  abline(h = 0.8, col = "red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  dev.off()
  
  r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  softPower <- if (any(r2 >= 0.8, na.rm = TRUE)) {
    sft$fitIndices[min(which(r2 >= 0.8)), 1]
  } else {
    sft$fitIndices[which.max(r2), 1]
  }
  
  net <- blockwiseModules(
    datExpr,
    power             = softPower,
    networkType       = "signed",
    TOMType           = "signed",
    minModuleSize     = minModuleSize,
    deepSplit         = deepSplit,
    reassignThreshold = 0,
    mergeCutHeight    = mergeCutHeight,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = paste0(prefix, "_TOM"),
    verbose           = 3
  )
  
  moduleLabels <- net$colors
  moduleColors <- labels2colors(net$colors)
  MEs          <- net$MEs
  save(MEs, moduleLabels, moduleColors, datExpr,
       file = paste0(prefix, "_network.RData"))
  
  module_df <- data.frame(Gene = colnames(datExpr),
                          Module = moduleLabels,
                          ModuleColor = moduleColors)
  module_df <- module_df[order(module_df$Module), ]
  write.csv(module_df, paste0(prefix, "_modules_genes.csv"), row.names = FALSE)
  
  genes_order <- module_df$Gene
  hm_mat <- datExpr[, genes_order, drop = FALSE]
  ann_col <- data.frame(
    Module = factor(module_df$Module, levels = sort(unique(module_df$Module))),
    row.names = genes_order
  )
  hm_mat_clipped <- hm_mat
  hm_mat_clipped[hm_mat_clipped >  3] <-  3
  hm_mat_clipped[hm_mat_clipped < -3] <- -3
  pdf(paste0(prefix, "_heatmap_modules_LFC_clip_-3to3.pdf"), width = 6, height = 4)
  pheatmap::pheatmap(hm_mat_clipped,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     scale = "none",
                     show_colnames = FALSE,
                     show_rownames = TRUE,
                     annotation_col = ann_col,
                     color = colorRampPalette(c("beige","orange"))(255),
                     main = "")
  dev.off()
  
  # Hallmark Gene Sets
  msig_h <- msigdbr(species = "Mus musculus", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  mod_list <- split(module_df$Gene, module_df$Module)
  mod_list <- mod_list[names(mod_list) != "0"]
  
  enrich_list <- lapply(names(mod_list), function(m) {
    g <- unique(mod_list[[m]])
    er <- tryCatch(
      enricher(g, TERM2GENE = msig_h,
               pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,
               qvalueCutoff = 0.05, pvalueCutoff = 1),
      error = function(e) NULL
    )
    if (!is.null(er)) {
      df <- as.data.frame(er)
      if (nrow(df) > 0) { df$Module <- m; return(df) }
    }
    NULL
  })
  enrich_df <- do.call(rbind, enrich_list)
  
  if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
    write.csv(enrich_df, paste0(prefix, "_hallmark_enrichment_ALL.csv"), row.names = FALSE)
    plot_df <- enrich_df |>
      transmute(
        Module   = factor(as.integer(Module), levels = sort(unique(as.integer(Module)))),
        Pathway  = ID,
        Count    = as.numeric(Count),
        padj     = as.numeric(p.adjust),
        GR_num   = as.numeric(sub("/.*","", GeneRatio)),
        GR_den   = as.numeric(sub(".*/","", GeneRatio)),
        GeneRatio= GR_num / GR_den,
        mlp      = -log10(padj)
      ) |>
      mutate(Pathway_pretty = sapply(Pathway, pretty_pathway_name))
    
    keep_paths <- plot_df |>
      group_by(Pathway_pretty) |>
      summarise(best = max(mlp, na.rm = TRUE), .groups = "drop") |>
      arrange(desc(best)) |>
      slice_head(n = 50) |>
      pull(Pathway_pretty)
    
    plot_df2 <- plot_df |>
      filter(Pathway_pretty %in% keep_paths) |>
      complete(Module, Pathway_pretty,
               fill = list(Count = 0, GeneRatio = 0, padj = NA_real_, mlp = NA_real_)) |>
      group_by(Pathway_pretty) |>
      mutate(path_order = max(mlp, na.rm = TRUE)) |>
      ungroup() |>
      mutate(Pathway_pretty = fct_reorder(Pathway_pretty, path_order))
    
    p <- ggplot(plot_df2, aes(x = Module, y = Pathway_pretty, size = GeneRatio, color = mlp)) +
      geom_point() +
      scale_size_continuous(range = c(1.2, 6), name = "GeneRatio") +
      scale_color_viridis_c(option = "C", na.value = "grey90", name = "-log10(adj. p)") +
      theme_bw(base_size = 10) +
      theme(axis.text.y = element_text(size = 7),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()) +
      labs(title = "", x = "WGCNA Module", y = "Pathway")
    
    ggsave(paste0(prefix, "_hallmark_dotplot_GLOBAL.pdf"), p, width = 10, height = 10)
    ggsave(paste0(prefix, "_hallmark_RC_dotplot_GLOBAL_q0.05.png"), p, width = 5, height = 4, dpi = 300)
  }
  
  invisible(list(datExpr = datExpr, module_df = module_df, net = net))
}

lfc_blck6 <- lfc_matrix_from_deglist(Blck6_up_Allo_vs_Syn, keep_ct = ct_keep)
lfc_balbc <- lfc_matrix_from_deglist(Balbc_up_Allo_vs_Syn, keep_ct = ct_keep)

res_blck6 <- run_wgcna_from_lfc(lfc_blck6, prefix = "Blck6_allo_LFC")
res_balbc <- run_wgcna_from_lfc(lfc_balbc, prefix = "Balbc_allo_LFC")

# Suppl. Fig. 13+14
ID_colors <- c(
  "BALB/c to BALB/c 1" = "grey70",
  "BALB/c to BALB/c 2" = "grey50",
  "C57BL/6 to C57BL/6 1" = "grey30",
  "BALB/c to C57BL/6 1" = "darkolivegreen2",
  "BALB/c to C57BL/6 2" = "darkolivegreen3",
  "C57BL/6 to BALB/c 1" = "lightblue3",
  "C57BL/6 to BALB/c 2" = "lightblue4",
  "C57BL/6 to BALB/c 3" = "lightblue"
)

ordered_ids <- c(
  "BALB/c to BALB/c 1","BALB/c to BALB/c 2","C57BL/6 to C57BL/6 1",
  "BALB/c to C57BL/6 1","BALB/c to C57BL/6 2",
  "C57BL/6 to BALB/c 1","C57BL/6 to BALB/c 2","C57BL/6 to BALB/c 3"
)

set_idents_plot <- function(obj, cols) {
  Idents(obj) <- "celltype_level_2"
  DimPlot(obj, label = FALSE, cols = cols) + NoLegend()
}

avg_heatmap <- function(obj, genes, order_vec) {
  avg <- AverageExpression(obj, features = rev(genes), group.by = "celltype_level_2", assays = "RNA", slot = "data")
  m <- avg$RNA
  m <- t(apply(m, 1, function(x) x / max(x, na.rm = TRUE)))
  colnames(m) <- gsub("-", "_", colnames(m))
  order_fixed <- gsub("-", "_", order_vec)
  m <- m[, order_fixed, drop = FALSE]
  pheatmap(
    m, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE,
    scale = "none", breaks = seq(0, 1, length.out = 101),
    color = colorRampPalette(c("black", "yellow"))(100),
    border_color = NA, legend = FALSE
  )
}

percent_by_id <- function(meta, id_col = "ID_to_plot", cell_col = "celltype_level_2") {
  meta %>%
    group_by(.data[[id_col]], .data[[cell_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[id_col]]) %>%
    mutate(total = sum(count), percent = count / total) %>%
    ungroup()
}

complete_grid <- function(df, id_col = "ID_to_plot", cell_col = "celltype_level_2", cell_levels = NULL) {
  ids <- unique(df[[id_col]])
  cells <- if (is.null(cell_levels)) unique(df[[cell_col]]) else cell_levels
  expand.grid(tmp_id = ids, tmp_cell = cells, stringsAsFactors = FALSE) %>%
    rename(!!id_col := tmp_id, !!cell_col := tmp_cell) %>%
    left_join(df, by = c(id_col, cell_col)) %>%
    mutate(count = ifelse(is.na(count), 0, count),
           total = ifelse(is.na(total), 0, total),
           percent = ifelse(is.na(percent), 0, percent))
}

bars_percent <- function(df, id_levels, cell_levels, id_colors, show_axes = TRUE) {
  df[[names(df)[1]]] <- factor(df[[names(df)[1]]], levels = id_levels)
  df[[names(df)[2]]] <- factor(df[[names(df)[2]]], levels = cell_levels)
  g <- ggplot(df, aes_string(x = names(df)[2], y = "percent", fill = names(df)[1])) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = id_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal()
  if (show_axes) {
    g + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    g + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
  }
}

bars_percent_groups <- function(df, group_col = "group", cell_col = "celltype_level_2") {
  df %>%
    group_by(.data[[group_col]], .data[[cell_col]]) %>%
    summarise(mean_percent = mean(percent), .groups = "drop") %>%
    ggplot(aes_string(x = cell_col, y = "mean_percent", fill = group_col)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
}

t_tests_by_celltype <- function(df, value_col = "percent", group_col = "group", cell_col = "celltype_level_2") {
  df %>%
    group_by(.data[[cell_col]]) %>%
    rstatix::t_test(reformulate(group_col, response = value_col)) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    select(!!cell_col, statistic, df, p, p.adj, group1, group2, n1, n2) %>%
    ungroup()
}

set_idents_plot(Ktx_DCT_mouse_subclustering, c("DCT1" = "grey40","DCT2" = "grey80","DCT_prolif" = "#36648B","DCT_Injury_m1" = "#FFB90F","DCT_Injury_m2" = "#EE2C2C"))
DCT_genes <- c("Slc12a3","Trpm7","Fnip2","Slc8a1","Trpm6","Calb1","Top2a","Mki67","Fth1","Ftl1","Gpx3","Mt1","H2-D1","Adamts1","Cd44","Samd5","Nfkbiz","Il34","Cxcl10","Osmr","Gbp3")
avg_heatmap(Ktx_DCT_mouse_subclustering, DCT_genes, c("DCT1","DCT2","DCT_prolif","DCT_Injury_m1","DCT_Injury_m2"))

df_mouse <- percent_by_id(Ktx_DCT_mouse_subclustering@meta.data, "ID_to_plot", "celltype_level_2")
df_mouse <- complete_grid(df_mouse, "ID_to_plot", "celltype_level_2", c("DCT1","DCT2","DCT_prolif","DCT_Injury_m1","DCT_Injury_m2"))
df_mouse$ID_to_plot <- factor(df_mouse$ID_to_plot, levels = ordered_ids)
df_mouse$celltype_level_2 <- factor(df_mouse$celltype_level_2, levels = c("DCT1","DCT2","DCT_prolif","DCT_Injury_m1","DCT_Injury_m2"))
bars_percent(df_mouse[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("DCT1","DCT2","DCT_prolif","DCT_Injury_m1","DCT_Injury_m2"), ID_colors, FALSE)

total_cells_by_ID_dct <- Ktx_DCT_mouse_subclustering@meta.data %>% group_by(ID) %>% summarise(total_cells = n(), .groups = "drop")
df_dct <- Ktx_DCT_mouse_subclustering@meta.data %>% group_by(ID, celltype_level_2) %>% summarise(count = n(), .groups = "drop") %>% left_join(total_cells_by_ID_dct, by = "ID") %>% mutate(percent = count / total_cells) %>% ungroup() %>% mutate(group = Ktx_DCT_mouse_subclustering$group[match(ID, Ktx_DCT_mouse_subclustering$ID)])
unique_IDs_dct <- df_dct %>% distinct(ID, group)
all_combinations_dct <- expand.grid(ID = unique_IDs_dct$ID, celltype_level_2 = unique(df_dct$celltype_level_2)) %>% left_join(unique_IDs_dct, by = "ID")
df_dct_complete <- all_combinations_dct %>% left_join(df_dct, by = c("ID","celltype_level_2","group")) %>% mutate(count = ifelse(is.na(count), 0, count), total_cells = ifelse(is.na(total_cells), 0, total_cells), percent = ifelse(is.na(percent), 0, percent))
t_test_results_dct <- t_tests_by_celltype(df_dct_complete, "percent", "group", "celltype_level_2")
print(t_test_results_dct)
df_dct_summary <- df_dct_complete %>% group_by(group, celltype_level_2) %>% summarise(mean_percent = mean(percent), .groups = "drop")
ggplot(df_dct_summary, aes(x = celltype_level_2, y = mean_percent, fill = group)) + geom_bar(stat = "identity", position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())

set_idents_plot(Ktx_CNT_mouse_subclustering, c("CNT1" = "darkseagreen","CNT2" = "darkseagreen1","CNT/CD" = "cornsilk3","CNT_prolif" = "#36648B","CNT_Injury_m1" = "#EE2C2C"))
CNT_genes <- c("Calb1","Trpm6","Slc8a1","Scnn1b","Scnn1g","Trpv5","S100g","Aqp2","Top2a","Mki67","Nfkbiz","Cxcl10","Gbp3","Samd5","Lcn2","Adamts1","Tpm1")
avg_heatmap(Ktx_CNT_mouse_subclustering, CNT_genes, c("CNT1","CNT2","CNT/CD","CNT_prolif","CNT_Injury_m1"))
df_cnt <- percent_by_id(Ktx_CNT_mouse_subclustering@meta.data, "ID_to_plot", "celltype_level_2") %>% complete_grid("ID_to_plot","celltype_level_2", c("CNT1","CNT2","CNT/CD","CNT_prolif","CNT_Injury_m1"))
df_cnt$ID_to_plot <- factor(df_cnt$ID_to_plot, levels = ordered_ids)
df_cnt$celltype_level_2 <- factor(df_cnt$celltype_level_2, levels = c("CNT1","CNT2","CNT/CD","CNT_prolif","CNT_Injury_m1"))
bars_percent(df_cnt[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("CNT1","CNT2","CNT/CD","CNT_prolif","CNT_Injury_m1"), ID_colors, FALSE)

set_idents_plot(Ktx_IntC_mouse_subclustering, c("Fibro1" = "#CDC8B1","Fibro2" = "#8B8878","Fibro3" = "#8EE5EE","Mixed" = "chocolate"))
IntC_genes <- c("Dapk2","Fbln5","Col14a1","Negr1","Spon1","Cdh11","Col1a1","Col1a2","Dcn","Vim","Cd44","Tnfrsf12a","Cxcl10","Cxcl9","Serpina3f","Rpl32","Top2a","Mki67","Gata3","Ren1","Acta2","Myh11")
avg_heatmap(Ktx_IntC_mouse_subclustering, IntC_genes, c("Fibro1","Fibro2","Fibro3","Mixed"))
df_intc <- percent_by_id(Integrated_IntC_mouse@meta.data, "ID_to_plot", "celltype_level_2") %>% complete_grid("ID_to_plot","celltype_level_2", c("Fibro1","Fibro2","Fibro3","REN/MC/VSMC"))
df_intc$ID_to_plot <- factor(df_intc$ID_to_plot, levels = ordered_ids)
df_intc$celltype_level_2 <- factor(df_intc$celltype_level_2, levels = c("Fibro1","Fibro2","Fibro3","REN/MC/VSMC"))
bars_percent(df_intc[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("Fibro1","Fibro2","Fibro3","REN/MC/VSMC"), ID_colors, TRUE)

set_idents_plot(Ktx_tL_mouse_subclustering, c("DTL1" = "orange","DTL2" = "orange3","DTL3" = "orangered","ATL1" = "#FFDEAD","ATL2" = "#CDB38B"))
tL_genes <- c("Cdh13","Slc9a3","Slc4a11","Slc14a2","Slit3","Cxcl10","Adamts1","Aqp1","Cacna1d","Sorcs3","Epha7","Kirrel3","Ephb2","Itga2","Clcnka","Cldn10","Sgcz")
avg_heatmap(Ktx_tL_mouse_subclustering, tL_genes, c("DTL1","DTL2","DTL3","ATL1","ATL2"))
df_tL <- percent_by_id(Integrated_tL_mouse@meta.data, "ID_to_plot", "celltype_level_2")
df_tL$ID_to_plot <- factor(df_tL$ID_to_plot, levels = ordered_ids)
df_tL$celltype_level_2 <- factor(df_tL$celltype_level_2, levels = c("DTL1","DTL2","DTL3","ATL1","ATL2"))
ggplot(df_tL, aes(x = celltype_level_2, y = percent, fill = ID_to_plot)) + geom_bar(stat = "identity", position = position_dodge()) + scale_fill_manual(values = ID_colors) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) + labs(y = "Relative Abundance (%)")
ggplot(df_tL, aes(x = celltype_level_2, y = percent, fill = ID_to_plot)) + geom_bar(stat = "identity", position = position_dodge()) + scale_fill_manual(values = ID_colors) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()

set_idents_plot(Ktx_EC_mouse_subclustering, c("Non-FenEC1" = "#FFC125","Non-FenEC2" = "orange3","Non-FenEC3" = "#EEE685","FenEC1" = "#8B864E","FenEC2" = "lemonchiffon3","EC_prolif" = "darkseagreen3","EC_Injury_m1" = "#EE6363","EC_Injury_m2" = "orangered"))
EC_genes <- c("Nebl","Fbln2","Slc14a1","Aqp1","Kdr","Plpp3","Plvap","Tll1","Plcl2","Top2a","Mki67","Diaph3","Ptprc","Arhgap15","Dock10","Ccr2","Serpine1","Tnfrsf23","Cxcl9","Cxcl10","Cd274","Ifi207","Gbp2","Iigp1","Vcam1","Tnfrsf22","Pecam1")
avg_heatmap(Ktx_EC_mouse_subclustering, EC_genes, c("Non-FenEC1","Non-FenEC2","Non-FenEC3","FenEC1","FenEC2","EC_prolif","EC_Injury_m1","EC_Injury_m2"))
df_ec <- percent_by_id(Ktx_EC_mouse_subclustering@meta.data, "ID_to_plot", "celltype_level_2") %>% complete_grid("ID_to_plot","celltype_level_2", c("Non-FenEC1","Non-FenEC2","Non-FenEC3","FenEC1","FenEC2","EC_prolif","EC_Injury_m1","EC_Injury_m2"))
df_ec$ID_to_plot <- factor(df_ec$ID_to_plot, levels = ordered_ids)
df_ec$celltype_level_2 <- factor(df_ec$celltype_level_2, levels = c("Non-FenEC1","Non-FenEC2","Non-FenEC3","FenEC1","FenEC2","EC_prolif","EC_Injury_m1","EC_Injury_m2"))
p1 <- bars_percent(df_ec[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("Non-FenEC1","Non-FenEC2","Non-FenEC3","FenEC1","FenEC2","EC_prolif","EC_Injury_m1","EC_Injury_m2"), ID_colors, TRUE)
p2 <- bars_percent(df_ec[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("Non-FenEC1","Non-FenEC2","Non-FenEC3","FenEC1","FenEC2","EC_prolif","EC_Injury_m1","EC_Injury_m2"), ID_colors, FALSE)
print(p1); print(p2)

set_idents_plot(Ktx_CD_PC_mouse_subclustering, c("cPC1" = "#FFF0F5","cPC2" = "lavenderblush3","mPC1" = "#8B3A62","mPC2" = "#CD6090","PC_prolif" = "lemonchiffon3"))
PC_genes <- c("Scnn1g","Scnn1b","Cd74","Cd274","Irf1","Pappa","Aqp3","Slc6a6","Ptgs1","Cdh16","Aqp2","Gata3","Slc14a2","Pcdh7","Fxyd2","Tacstd2","Mif","Scarb1","Pdgfb","Igf1r","Tnfrsf19","Nfkbia","Tnfaip2","Pfkb3","Top2a","Mki67","Brca1")
avg_heatmap(Ktx_CD_PC_mouse_subclustering, PC_genes, c("cPC1","cPC2","mPC1","mPC2","PC_prolif"))
df_pc <- percent_by_id(Integrated_PC_mouse@meta.data, "ID_to_plot", "celltype_level_2") %>% complete_grid("ID_to_plot","celltype_level_2", c("cPC1","cPC2","mPC1","mPC2","PC_prolif"))
df_pc$ID_to_plot <- factor(df_pc$ID_to_plot, levels = ordered_ids)
df_pc$celltype_level_2 <- factor(df_pc$celltype_level_2, levels = c("cPC1","cPC2","mPC1","mPC2","PC_prolif"))
bars_percent(df_pc[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("cPC1","cPC2","mPC1","mPC2","PC_prolif"), ID_colors, FALSE)

set_idents_plot(Ktx_CD_IC_mouse_subclustering, c("CD-IC-A1" = "#0000FF","CD-IC-A2" = "#00008B","CD-IC-B1" = "#5F9EA0","CD-IC-B2" = "#98F5FF"))
IC_genes <- c("Slc26a7","Aqp6","Scnn1b","Slc4a1","Ptprd","Oxgr1","Slc26a4","Slc4a9","Tmem72","Slc9a9","Casr","Mif","Lcp2")
avg_heatmap(Ktx_CD_IC_mouse_subclustering, IC_genes, c("CD-IC-A1","CD-IC-A2","CD-IC-B1","CD-IC-B2"))
df_ic <- percent_by_id(Ktx_CD_IC_mouse_subclustering@meta.data, "ID_to_plot", "celltype_level_2") %>% complete_grid("ID_to_plot","celltype_level_2", c("CD-IC-A1","CD-IC-A2","CD-IC-B1","CD-IC-B2"))
df_ic$ID_to_plot <- factor(df_ic$ID_to_plot, levels = ordered_ids)
df_ic$celltype_level_2 <- factor(df_ic$celltype_level_2, levels = c("CD-IC-A1","CD-IC-A2","CD-IC-B1","CD-IC-B2"))
bars_percent(df_ic[,c("ID_to_plot","celltype_level_2","percent")], ordered_ids, c("CD-IC-A1","CD-IC-A2","CD-IC-B1","CD-IC-B2"), ID_colors, FALSE)


# Suppl. Fig. 19
Ktx_data_human <- readRDS(".../Ktx_data_human.rds")
Idents(Ktx_data_human) <- "sample"

ordered_ids <- c("C1", "C2", "C4", "R1", "R3", "R4")
Ktx_data_human$sample <- factor(Ktx_data_human$sample, levels = ordered_ids)

sample_colors <- c("C1" = "white", "C2" = "white", "C4" = "white", 
                   "R1" = "firebrick3", "R3" = "firebrick3", "R4" = "firebrick3")

features_with_limits <- list(
  "nCount_RNA" = 60000,  
  "nFeature_RNA" = 12000, 
  "percent.mt" = 5     
)

for (feature in names(features_with_limits)) {
  feature_y_limit <- features_with_limits[[feature]]
  
  p <- VlnPlot(
    Ktx_data_human, 
    features = c(feature), 
    pt.size = 0,
    group.by = "sample"  
  ) + 
    scale_fill_manual(values = sample_colors) +
    scale_x_discrete(limits = ordered_ids) + 
    theme(
      axis.text.x = element_text(size = rel(0.8))
    ) +
    NoLegend()
  
  if (!is.null(feature_y_limit)) {
    p <- p + coord_cartesian(ylim = c(0, feature_y_limit))
  }
  
  print(p)
}


cell_counts <- Ktx_data_human@meta.data %>%
  group_by(sample) %>%
  summarise(cell_count = n())

ggplot(cell_counts, aes(x = sample, y = cell_count)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 12),  
    axis.text.y = element_text(size = 10)
  ) +
  ylab("")
