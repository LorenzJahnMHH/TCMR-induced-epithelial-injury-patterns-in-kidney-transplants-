# Figure 2

library(Seurat)
library(tidyverse)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

setwd("D:/R_Ordner/KTx Daten")
Ktx_data <- readRDS("Ktx_data.rds")

setwd("D:/R_Ordner/Ktx Daten/Blck6_Balbc_DGE_lists")

cell_types <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "IntC", "Uro", "PEC", "Leuko")

Balbc_Allo_marker <- list()
Balbc_up_Allo_vs_Syn <- list()
Balbc_down_Allo_vs_Syn <- list()

Blck6_Allo_marker <- list()
Blck6_up_Allo_vs_Syn <- list()
Blck6_down_Allo_vs_Syn <- list()

for (cell_type in cell_types) {
  
  ktx_cells <- subset(Ktx_data, subset = celltype_level_1 == cell_type)
  
  DefaultAssay(ktx_cells) <- "RNA"
  Idents(ktx_cells) <- "groupID"
  
  Balbc_Allo_marker[[cell_type]] <- FindMarkers(ktx_cells, ident.1 = 'Balbc_Blck6_alo', ident.2 = 'Balbc_Balbc_syn', logfc.threshold = 1)
  

  Balbc_up_Allo_vs_Syn[[cell_type]] <- Balbc_Allo_marker[[cell_type]] %>%
    filter(p_val_adj < 0.05, avg_log2FC > 1, pct.1 > 0.1)
  

  Balbc_down_Allo_vs_Syn[[cell_type]] <- Balbc_Allo_marker[[cell_type]] %>%
    filter(p_val_adj < 0.05, avg_log2FC < -1, pct.2 > 0.1)
  
  Blck6_Allo_marker[[cell_type]] <- FindMarkers(ktx_cells, ident.1 = 'Blck6_Balbc_alo', ident.2 = 'Blck6_Blck6_syn', logfc.threshold = 1)
  
  Blck6_up_Allo_vs_Syn[[cell_type]] <- Blck6_Allo_marker[[cell_type]] %>%
    filter(p_val_adj < 0.05, avg_log2FC > 1, pct.1 > 0.1)
  
  Blck6_down_Allo_vs_Syn[[cell_type]] <- Blck6_Allo_marker[[cell_type]] %>%
    filter(p_val_adj < 0.05, avg_log2FC < -1, pct.2 > 0.1)
  
  rm(ktx_cells)
  gc()
  

  if (nrow(Balbc_up_Allo_vs_Syn[[cell_type]]) > 0) {
    Balbc_up_Allo_vs_Syn[[cell_type]]$cell_type <- paste0("Balbc_", cell_type)
    Balbc_up_Allo_vs_Syn[[cell_type]]$genes <- rownames(Balbc_up_Allo_vs_Syn[[cell_type]])
  }
  
  if (nrow(Balbc_down_Allo_vs_Syn[[cell_type]]) > 0) {
    Balbc_down_Allo_vs_Syn[[cell_type]]$cell_type <- paste0("Balbc_", cell_type)
    Balbc_down_Allo_vs_Syn[[cell_type]]$genes <- rownames(Balbc_down_Allo_vs_Syn[[cell_type]])
  }
  
  if (nrow(Blck6_up_Allo_vs_Syn[[cell_type]]) > 0) {
    Blck6_up_Allo_vs_Syn[[cell_type]]$cell_type <- paste0("Blck6_", cell_type)
    Blck6_up_Allo_vs_Syn[[cell_type]]$genes <- rownames(Blck6_up_Allo_vs_Syn[[cell_type]])
  }
  
  if (nrow(Blck6_down_Allo_vs_Syn[[cell_type]]) > 0) {
    Blck6_down_Allo_vs_Syn[[cell_type]]$cell_type <- paste0("Blck6_", cell_type)
    Blck6_down_Allo_vs_Syn[[cell_type]]$genes <- rownames(Blck6_down_Allo_vs_Syn[[cell_type]])
  }
}


saveRDS(Balbc_up_Allo_vs_Syn, file = "Balbc_up_Allo_vs_Syn.rds")
saveRDS(Balbc_down_Allo_vs_Syn, file = "Balbc_down_Allo_vs_Syn.rds")
saveRDS(Blck6_up_Allo_vs_Syn, file = "Blck6_up_Allo_vs_Syn.rds")
saveRDS(Blck6_down_Allo_vs_Syn, file = "Blck6_down_Allo_vs_Syn.rds")


save_to_excel <- function(data_list, filename) {
  wb <- createWorkbook()
  
  for (cell_type in names(data_list)) {
    addWorksheet(wb, sheetName = cell_type)
    writeData(wb, sheet = cell_type, data_list[[cell_type]], rowNames = TRUE)
  }
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}


save_to_excel(Balbc_up_Allo_vs_Syn, "Balbc_up_Allo_vs_Syn.xlsx")
save_to_excel(Balbc_down_Allo_vs_Syn, "Balbc_down_Allo_vs_Syn.xlsx")
save_to_excel(Blck6_up_Allo_vs_Syn, "Blck6_up_Allo_vs_Syn.xlsx")
save_to_excel(Blck6_down_Allo_vs_Syn, "Blck6_down_Allo_vs_Syn.xlsx")


# Figure 2A
Balbc_up_Allo_vs_Syn <- readRDS("Balbc_up_Allo_vs_Syn.rds")
Balbc_down_Allo_vs_Syn <- readRDS("Balbc_down_Allo_vs_Syn.rds")
Blck6_up_Allo_vs_Syn <- readRDS("Blck6_up_Allo_vs_Syn.rds")
Blck6_down_Allo_vs_Syn <- readRDS("Blck6_down_Allo_vs_Syn.rds")

count_genes <- function(data_list) {
  sapply(data_list, nrow)  
}

upregulated_counts_Blck6 <- count_genes(Blck6_up_Allo_vs_Syn)
downregulated_counts_Blck6 <- count_genes(Blck6_down_Allo_vs_Syn)
upregulated_counts_Balbc <- count_genes(Balbc_up_Allo_vs_Syn)
downregulated_counts_Balbc <- count_genes(Balbc_down_Allo_vs_Syn)

# Leukocytes are excluded due to insufficient cell counts in syngeneic samples for reliable DGE analysis
cell_types <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "IntC", "Uro", "PEC")

upregulated_counts_Blck6 <- upregulated_counts_Blck6[cell_types]
downregulated_counts_Blck6 <- downregulated_counts_Blck6[cell_types]
upregulated_counts_Balbc <- upregulated_counts_Balbc[cell_types]
downregulated_counts_Balbc <- downregulated_counts_Balbc[cell_types]

data <- data.frame(
  Cell_Type = rep(cell_types, times = 4),
  Count = c(upregulated_counts_Blck6, downregulated_counts_Blck6,
            upregulated_counts_Balbc, downregulated_counts_Balbc),
  Strain = rep(c("Blck6", "Blck6", "Balbc", "Balbc"), each = length(cell_types)),
  Regulation = rep(c("Upregulated", "Downregulated", "Upregulated", "Downregulated"), each = length(cell_types))
)

data$Interaction <- interaction(data$Strain, data$Regulation)
data$Interaction <- factor(data$Interaction, levels = c("Blck6.Upregulated", "Blck6.Downregulated", 
                                                        "Balbc.Upregulated", "Balbc.Downregulated"))

data$Cell_Type <- factor(data$Cell_Type, levels = cell_types)

upregulated_data <- subset(data, Regulation == "Upregulated")
downregulated_data <- subset(data, Regulation == "Downregulated")

ggplot(upregulated_data, aes(x = Cell_Type, y = Count, fill = Strain)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) + 
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),  
    axis.text.y = element_text(size = 13), 
    axis.title.x = element_text(size = 15),  
    axis.title.y = element_text(size = 15),  
    legend.position = "right"        
  ) + 
  ylim(0, max(upregulated_data$Count) * 1.1)

ggplot(downregulated_data, aes(x = Cell_Type, y = Count, fill = Strain)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) + 
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),  
    axis.text.y = element_text(size = 13), 
    axis.title.x = element_text(size = 15),  
    axis.title.y = element_text(size = 15),  
    legend.position = "right"        
  ) + 
  ylim(0, max(downregulated_data$Count) * 1.1)

ggplot(upregulated_data, aes(x = Cell_Type, y = Count, fill = Strain)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) + 
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(), 
    legend.position = "none",       
    panel.grid.minor = element_blank()  
  ) + 
  ylim(0, max(upregulated_data$Count) * 1.1)

ggplot(downregulated_data, aes(x = Cell_Type, y = Count, fill = Strain)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) + 
  scale_fill_manual(values = c("Blck6" = "skyblue", "Balbc" = "lightgreen")) + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(), 
    legend.position = "none",        
    panel.grid.minor = element_blank() 
  ) + 
  ylim(0, max(downregulated_data$Count) * 1.1)


# Figure 2B
setwd("D:/R_Ordner/Ktx Daten/Blck6_Balbc_DGE_lists")
# .../Blck6_Balbc_DGE_lists/


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
  show_rownames = TRUE, # Zeige die Gen-Namen
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
  show_rownames = TRUE, # Zeige die Gen-Namen
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
  show_rownames = FALSE, # Verstecke die Gen-Namen
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
  show_rownames = FALSE, # Verstecke die Gen-Namen
  show_colnames = FALSE,
  labels_col = NULL,
  legend = FALSE
)

# Fig 2C
# "_pathways_raw.txt"- files were generated by computing DGE gene overlaps with Hallmark database on https://maayanlab.cloud/Enrichr/ with the DGE genes showed in Fig 2A
#example name of a resulting txt. file: "Blck6_PT_pathways_raw.txt"

balbc_dir <- "D:/R_Ordner/KTx Daten/Pathway_enrichment/Balbc_up_pct_01/"
blck6_dir <- "D:/R_Ordner/KTx Daten/Pathway_enrichment/Blck6_up_pct_01/"

process_enrichr_file <- function(file_path) {
  data <- read.delim(file_path, header = TRUE, sep = "\t", fill = TRUE, quote = "", check.names = FALSE)
  clean_data <- data %>%
    dplyr::select(Term, `Adjusted P-value`) %>%
    dplyr::rename(Pathway = Term, AdjustedPValue = `Adjusted P-value`) %>%
    dplyr::mutate(AdjustedPValue = as.numeric(AdjustedPValue))
  return(clean_data)
}

celltypes <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "IntC", "Uro", "PEC")

all_pathway_data <- data.frame(Pathway = character(),
                               AdjustedPValue = numeric(),
                               CellType = character(),
                               Strain = character(),
                               stringsAsFactors = FALSE)

for (ct in celltypes) {
  balbc_file <- paste0(balbc_dir, "Balbc_", ct, "_pathways_raw.txt")
  blck6_file <- paste0(blck6_dir, "Blck6_", ct, "_pathways_raw.txt")
  
  if (file.exists(balbc_file)) {
    balbc_data <- process_enrichr_file(balbc_file)
    balbc_data$CellType <- ct
    balbc_data$Strain <- "Balbc"
    all_pathway_data <- rbind(all_pathway_data, balbc_data)
  }
  
  if (file.exists(blck6_file)) {
    blck6_data <- process_enrichr_file(blck6_file)
    blck6_data$CellType <- ct
    blck6_data$Strain <- "Blck6"
    all_pathway_data <- rbind(all_pathway_data, blck6_data)
  }
}

filtered_data <- all_pathway_data %>%
  dplyr::filter(AdjustedPValue < 0.05) %>%
  mutate(log10_fdr = -log10(AdjustedPValue))

reshaped_data <- filtered_data %>%
  unite("CellType_Strain", CellType, Strain, sep = "_") %>%
  dcast(Pathway ~ CellType_Strain, value.var = "log10_fdr", fill = 0)

heatmap_matrix <- as.matrix(reshaped_data[, -1])
rownames(heatmap_matrix) <- reshaped_data$Pathway

pathway_counts <- filtered_data %>%
  group_by(Pathway) %>%
  summarize(celltype_count = sum(log10_fdr > 0, na.rm = TRUE)) %>%
  arrange(celltype_count)

top_pathways <- filtered_data %>%
  group_by(Pathway) %>%
  summarize(min_adj_pval = min(AdjustedPValue, na.rm = TRUE)) %>%
  arrange(min_adj_pval) %>%
  slice_head(n = 25)

filtered_top_data <- filtered_data %>%
  filter(Pathway %in% top_pathways$Pathway)

filtered_top_data$Pathway <- factor(filtered_top_data$Pathway, levels = pathway_counts$Pathway)

filtered_top_data$CellType <- factor(filtered_top_data$CellType, levels = celltypes)

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

