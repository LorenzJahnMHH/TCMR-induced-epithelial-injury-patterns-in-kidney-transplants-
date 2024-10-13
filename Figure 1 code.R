# Figure 1

setwd("D:/R_Ordner/Ktx Daten")
Ktx_data <- readRDS("Ktx_data.rds")

# necessary file : "Ktx_data.rds"

library(Seurat)
library(ggplot2)
library(dplyr)



# Figure 1C UMAP
Idents(Ktx_data) <- "celltype_level_1"

DimPlot(Ktx_data, label = F, cols = c('PT' = '#FF7256', 'tL' = '#8B2323', 'TAL' = '#00008B', 'PEC' = '#FAEBD7', 
                                      'Leuko' = '#9370DB', 'Podo' = '#FF7F00', 
                                      'IntC' = '#DDA0DD', 'EC' = '#B0E2FF', 'CD-PC' = '#9ACD32', 
                                      'CNT' = '#6B8E23', 'DCT' = '#BDB76B', 'Prolif' = '#00C5CD', 'Uro' = '#FFC125', 
                                      'CD-IC-A' = '#CDCDC1', 'CD-IC-B' = '#8B8B83')) + NoLegend() 

#Figure 1C Heatmap
Idents(Ktx_data) <- "celltype_level_1"

top_markers <- FindAllMarkers(
  object = Ktx_data, 
  only.pos = TRUE, 
  min.pct = 0.25,
  logfc.threshold = 0.5  
)

celltype_order = c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", 
                   "CD-IC-A", "CD-IC-B", "EC", "Leuko", "IntC", 
                   "Uro", "PEC", "Prolif")

top50 <- top_markers %>%
  filter(pct.1 > 0.5) %>%  
  group_by(cluster) %>%
  arrange(p_val_adj) %>%  
  slice_head(n = 50) %>%  
  ungroup()


top_genes <- top50$gene
top_genes_order <- top50 %>% 
  mutate(cluster = factor(cluster, levels = rev(celltype_order))) %>%  
  arrange(cluster) %>% 
  pull(gene)

top_genes_order <- rev(top_genes_order)

cpm_matrix <- matrix(0, nrow = length(top_genes_order), ncol = length(celltype_order))
rownames(cpm_matrix) <- top_genes_order
colnames(cpm_matrix) <- celltype_order


for (ct in celltype_order) {
  cells.tmp <- rownames(Ktx_data@meta.data[Ktx_data$celltype_level_1 == ct, ])
  
  if (length(cells.tmp) == 0) {
    next  
  }
  
  gene_counts <- Ktx_data@assays$RNA@counts[top_genes_order, cells.tmp, drop = FALSE]
  gene_counts <- as.matrix(gene_counts)
  total_counts <- colSums(gene_counts)
  total_counts[total_counts == 0] <- 1
  
  cpm <- t(t(gene_counts) / total_counts * 1e6)
  cpm_matrix[, ct] <- rowMeans(cpm, na.rm = TRUE)
}


cpm_matrix_norm <- t(apply(cpm_matrix, 1, function(x) x / max(x, na.rm = TRUE)))

cpm_df <- as.data.frame(as.table(cpm_matrix_norm))
colnames(cpm_df) <- c("Gene", "CellType", "CPM")

cpm_df$Gene <- factor(cpm_df$Gene, levels = unique(top_genes_order))
cpm_df$CellType <- factor(cpm_df$CellType, levels = celltype_order)

p <- ggplot(cpm_df, aes(x = CellType, y = Gene, fill = CPM)) + 
  geom_tile(color = NA) +  
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow", 
                       na.value = "black", limits = c(0, 1)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),  
    axis.ticks = element_blank(),  
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1)
  )

print(p)


p <- ggplot(cpm_df, aes(x = CellType, y = Gene, fill = CPM)) + 
  geom_tile(color = NA) +  
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow", 
                       na.value = "black", limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(),  
    panel.grid = element_blank(),   
    panel.border = element_blank(),
    legend.position = "none"  
  ) 

print(p)


# Figure 1E Barplot
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

df_mouse <- Ktx_data@meta.data %>%
  group_by(ID_to_plot, celltype_level_1) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(ID_to_plot) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_ids <- c("BALB/c to BALB/c 1", "BALB/c to BALB/c 2", "C57BL/6 to C57BL/6 1",
                 "BALB/c to C57BL/6 1", "BALB/c to C57BL/6 2", 
                 "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3")

df_mouse$ID_to_plot <- factor(df_mouse$ID_to_plot, levels = ordered_ids)

celltype_order <- c("Podo", "PT", "tL", "TAL", "MD", "DCT", "CNT", "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "Leuko", "IntC", "Uro", "PEC", "Prolif")

df_mouse$celltype_level_1 <- factor(df_mouse$celltype_level_1, levels = celltype_order)

df_mouse_long <- df_mouse %>%
  select(ID_to_plot, celltype_level_1, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

ID_colors <- c("BALB/c to BALB/c 1" = "grey70", 
               "BALB/c to BALB/c 2" = "grey50", 
               "C57BL/6 to C57BL/6 1" = "grey30",
               "BALB/c to C57BL/6 1" = "darkolivegreen2",
               "BALB/c to C57BL/6 2" = "darkolivegreen3",
               "C57BL/6 to BALB/c 1" = "lightblue3",
               "C57BL/6 to BALB/c 2" = "lightblue4",
               "C57BL/6 to BALB/c 3" = "lightblue")

ggplot(df_mouse_long, aes(x = celltype_level_1, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cell Type", y = "Relative Abundance (%)", fill = "Sample") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 13), 
    legend.position = "right" 
  )

 ggplot(df_mouse_long, aes(x = celltype_level_1, y = value, fill = ID_to_plot)) +
geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 13), 
    axis.text.y.left  = element_blank(),
    legend.position = "none" 
  )

