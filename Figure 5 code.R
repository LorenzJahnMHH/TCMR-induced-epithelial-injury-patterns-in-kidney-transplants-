# Figure 5

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(biomaRt)

# Figure 5A UMAP
Ktx_data_human <- readRDS("Ktx_human_int_harmony.rds")

Idents(Ktx_data_human) <- "celltype_level_1"
DimPlot(Ktx_data_human, label = F, cols = c('PT' = 'coral1', 'tL' = 'brown4', 'TAL' = 'darkblue', 'PEC' = 'antiquewhite', 
                                                   'Leuko' = 'mediumpurple', 'Podo' = 'darkorange1', 
                                                   'IntC' = 'plum', 'EC' = 'lightskyblue1', 'CD-PC' = 'olivedrab3', 
                                                   'CNT' = 'olivedrab', 'DCT' = 'darkkhaki', 'Uro' = 'goldenrod1', 
                                                   'CD-IC' = 'ivory3', 'Myo' = 'ivory4'))    + NoLegend()

# Figure 5A Heatmap
Idents(Ktx_data_human) <- "celltype_level_1"

top_markers <- FindAllMarkers(
  object = Ktx_data_human, 
  only.pos = TRUE, 
  min.pct = 0.5,
  logfc.threshold = 0.5  
)

celltype_order <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", 
                    "CD-IC", "EC", "Leuko", "IntC", 
                    "Uro", "PEC", "Myo")

top50 <- top_markers %>%
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

avg_expression <- AverageExpression(
  Ktx_data_human,
  features = top_genes_order,
  group.by = "celltype_level_1",
  assays = "RNA"
)$RNA

avg_expression <- avg_expression[top_genes_order, celltype_order, drop = FALSE]
normalized_matrix <- t(apply(avg_expression, 1, function(x) x / max(x, na.rm = TRUE)))

heatmap_df <- as.data.frame(as.table(normalized_matrix))
colnames(heatmap_df) <- c("Gene", "CellType", "Expression")

heatmap_df$Gene <- factor(heatmap_df$Gene, levels = unique(top_genes_order))
heatmap_df$CellType <- factor(heatmap_df$CellType, levels = celltype_order)


p <- ggplot(heatmap_df, aes(x = CellType, y = Gene, fill = Expression)) + 
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


#Figure 5B
df_human <- Ktx_data_human@meta.data %>%
  group_by(sample, celltype_level_1) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_samples <- c("C1", "C2", "C4", "R1", "R3", "R4")
df_human$sample <- factor(df_human$sample, levels = ordered_samples)

celltype_order <- c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", 
                    "CD-IC", "EC", "Leuko", "IntC", "Uro", "PEC", "Myo")
df_human$celltype_level_1 <- factor(df_human$celltype_level_1, levels = celltype_order)

sample_colors <- c("C1" = "grey70", 
                   "C2" = "grey50", 
                   "C4" = "grey30",
                   "R1" = "lightcoral",
                   "R3" = "indianred3",
                   "R4" = "brown3")

p1 <- ggplot(df_human, aes(x = celltype_level_1, y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cell Type", y = "Relative Abundance (%)", fill = "Sample") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 13), 
    legend.position = "right" 
  )

p2 <- ggplot(df_human, aes(x = celltype_level_1, y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
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

print(p1)

print(p2)



# Figure 5C UMAP
PT_human_KTX_subclustering <- readRDS("PT_human_KTX_subclustering.rds")
TAL_human_KTX_subclustering <- readRDS("TAL_human_KTX_subclustering.rds")

DimPlot(PT_human_KTX_subclustering, group.by = 'celltype_level_2', label = F, 
        cols = c('PT_S1' = 'coral3',
                 'PT_S2' = 'firebrick',  
                 'PT_S3' = 'coral4', 
                 'PT_S3_medullary' = 'darkorange2',
                 'PT_S1/2' = 'grey25', 
                 'PT_Injury_h1' = 'darkgoldenrod1', 
                 'PT_Injury_h2' = 'firebrick2'))  + NoLegend()

#Figure 5C Heatmap
desired_order <- c("PT S1", "PT S1/2", "PT S2",  "PT S3", "PT S3 med.", "PT Injury h1", "PT Injury h2")

gene_order <- rev(c("SLC5A2", "SLC7A7", "SLC22A8", "SLC34A1", "SLC7A13", "SLC22A24", "SLC22A7", "SLC25A47",
                    "UQCRQ", "ATP5F1E", "FTL", "FTH1", "SERPINA1",  
                   "RPL5", "RPL30", "SPARC", "HLA-A", "IFITM2", "KLF6", "VIM", "CD44", "VCAM1"))

avg_exp <- AverageExpression(
  PT_human_KTX_subclustering, 
  features = gene_order, 
  group.by = "celltype_to_plot", 
  assays = "RNA", 
  slot = "data"  
)

avg_exp_matrix <- avg_exp$RNA[gene_order, desired_order]

avg_exp_matrix_norm <- t(apply(avg_exp_matrix, 1, function(x) x / max(x, na.rm = TRUE)))

pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)

#Figure 5D
df_human <- PT_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_samples <- c("C1", "C2", "C4", "R1", "R3", "R4")

df_human$sample <- factor(df_human$sample, levels = ordered_samples)

celltype_order <- c("PT_S1", "PT_S1/2", "PT_S2", "PT_S3", "PT_S3_medullary", "PT_Injury_h1", "PT_Injury_h2")

df_human$celltype_level_2 <- factor(df_human$celltype_level_2, levels = celltype_order)

df_human_long <- df_human %>%
  select(sample, celltype_level_2, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

sample_colors <- c("C1" = "grey70", 
                   "C2" = "grey50", 
                   "C4" = "grey30",
                   "R1" = "lightcoral",
                   "R3" = "indianred",
                   "R4" = "firebrick")

ggplot(df_human_long, aes(x = celltype_level_2, y = value, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 13) 
  )

ggplot(df_human_long, aes(x = celltype_level_2, y = value, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),    
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),   
    legend.position = "none"         
  )



#Figure 5E 
Idents(PT_human_KTX_subclustering) <- "celltype_to_plot"  
Idents(Ktx_PT_mouse_subclustering) <- "celltype_to_plot"  

PT_mouse_marker <- FindAllMarkers(
  Ktx_PT_mouse_subclustering, 
  only.pos = TRUE, 
  logfc.threshold = 1
)
PT_mouse_marker <- PT_mouse_marker %>%
  filter(p_val_adj < 0.05, pct.1 > 0.1)

all_PT_genes <- unique(PT_mouse_marker$gene)

convertMouseGeneList <- function(x) {
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x,
                    mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = F)
  return(genesV2)
}

gene_mapping <- convertMouseGeneList(all_PT_genes)

gene_mapping <- gene_mapping[!duplicated(gene_mapping$MGI.symbol),]
gene_mapping <- gene_mapping[!duplicated(gene_mapping$HGNC.symbol),]

mouse_genes_to_scale <- gene_mapping$MGI.symbol
human_genes_to_scale <- gene_mapping$HGNC.symbol

Ktx_PT_mouse_subclustering <- ScaleData(Ktx_PT_mouse_subclustering, features = mouse_genes_to_scale)
PT_human_KTX_subclustering <- ScaleData(PT_human_KTX_subclustering, features = human_genes_to_scale)

mouse_agg <- AverageExpression(
  Ktx_PT_mouse_subclustering, 
  assays = "RNA", 
  features = mouse_genes_to_scale, 
  group.by = "celltype_to_plot", 
  slot = "scale.data"
)

human_agg <- AverageExpression(
  PT_human_KTX_subclustering, 
  assays = "RNA", 
  features = human_genes_to_scale, 
  group.by = "celltype_to_plot",  
  slot = "scale.data"
)

gene_mapping_filtered <- gene_mapping[gene_mapping$MGI.symbol %in% rownames(mouse_agg[["RNA"]]) & 
                                        gene_mapping$HGNC.symbol %in% rownames(human_agg[["RNA"]]), ]

mouse_exp <- as.matrix(mouse_agg[["RNA"]][gene_mapping_filtered$MGI.symbol, ])
human_exp <- as.matrix(human_agg[["RNA"]][gene_mapping_filtered$HGNC.symbol, ])

cor_matrix <- cor(mouse_exp, human_exp, method = "spearman")


breaks = seq(0, 0.5, by = 0.02)

pheatmap(
  cor_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = F, 
  color = colorRampPalette(c("white", "firebrick3"))(length(breaks) - 1),  
  breaks = breaks,  
  legend = TRUE, 
  fontsize = 12,  
  fontsize_row = 12,  
  fontsize_col = 12   
)


#Figure 5F UMAP
Idents(TAL_human_KTX_subclustering ) <- "celltype_level_2"
DimPlot(TAL_human_KTX_subclustering , group.by = 'celltype_level_2', label = F, 
        cols = c('cTAL1' = 'royalblue', 
                 'cTAL2' = 'royalblue4',
                 'mTAL1' = 'chartreuse3', 
                 'mTAL2' = 'darkseagreen',
                 'mTAL3' = 'darkgreen', 
                 'cTAL3' = 'lightskyblue2', 
                 'cTAL4' = 'lightskyblue4', 
                 'TAL_Injury_h1' = 'mediumpurple1',
                 'TAL_Injury_h2' = 'firebrick3')) + NoLegend()  

#Figure 5F Heatmap
desired_order_tal <- c("cTAL1", "cTAL2", "cTAL3", "cTAL4", "mTAL1", "mTAL2", "mTAL3", "TAL Injury h1", "TAL Injury h2")

gene_order_tal <- rev(c(
  
  "CLDN16", "TMEM52B", "KCNJ10",
  "PHACTR1", "CABP1",
  "CALCR", "FGL2", "PDE3A",  
  "RDH8", "PCBP3", "FMN1",
  "SLC22A2", "SLC25A5", "ATP1A1",
  "KLF2", "HLA-A", "S100A6", "MET", 
  "ITGB1", "ERO1A", "SLC2A1",  "TPM1", "SPP1", "ADAMTS1", "NFKBIZ", "CD44"
  
))

avg_exp_tal <- AverageExpression(
  TAL_human_KTX_subclustering, 
  features = gene_order_tal, 
  group.by = "celltype_to_plot", 
  assays = "RNA", 
  slot = "data"  
)

avg_exp_matrix_tal <- avg_exp_tal$RNA[gene_order_tal, desired_order_tal]

avg_exp_matrix_tal_norm <- t(apply(avg_exp_matrix_tal, 1, function(x) x / max(x, na.rm = TRUE)))

pheatmap(
  avg_exp_matrix_tal_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)


#Figure 5G
library(Seurat)
library(dplyr)
library(ggplot2)

desired_order <- c('cTAL1', 'cTAL2', 'cTAL3', 'cTAL4', 
                   'mTAL1', 'mTAL2', 'mTAL3', 
                   'TAL_Injury_h1', 'TAL_Injury_h2')

df_tal <- TAL_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

df_tal$celltype_level_2 <- factor(df_tal$celltype_level_2, levels = desired_order)

sample_colors <- c("C1" = "grey70", 
                   "C2" = "grey50", 
                   "C4" = "grey30",
                   "R1" = "lightcoral",
                   "R3" = "indianred",
                   "R4" = "firebrick")

ggplot(df_tal, aes(x = celltype_level_2, y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

ggplot(df_tal, aes(x = celltype_level_2, y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),    
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),   
    legend.position = "none"          
  )



#Figure 5H Correlation heatmap
Idents(TAL_human_KTX_subclustering) <- "celltype_to_plot"  
Idents(Ktx_TAL_mouse_subclustering) <- "celltype_to_plot" 

TAL_marker <- FindAllMarkers(
  Ktx_TAL_mouse_subclustering, 
  only.pos = TRUE, 
  logfc.threshold = 1
)
TAL_marker <- TAL_marker %>%
  filter(p_val_adj < 0.05, pct.1 > 0.1)

all_TAL_genes <- unique(TAL_marker$gene)

gene_mapping_TAL <- convertMouseGeneList(all_TAL_genes)

gene_mapping_TAL <- gene_mapping_TAL[!duplicated(gene_mapping_TAL$MGI.symbol),]
gene_mapping_TAL <- gene_mapping_TAL[!duplicated(gene_mapping_TAL$HGNC.symbol),]

mouse_genes_to_scale <- gene_mapping_TAL$MGI.symbol
human_genes_to_scale <- gene_mapping_TAL$HGNC.symbol

Ktx_TAL_mouse_subclustering <- ScaleData(Ktx_TAL_mouse_subclustering, features = mouse_genes_to_scale)
TAL_human_KTX_subclustering <- ScaleData(TAL_human_KTX_subclustering, features = human_genes_to_scale)

mouse_agg_TAL <- AverageExpression(
  Ktx_TAL_mouse_subclustering, 
  assays = "RNA", 
  features = mouse_genes_to_scale, 
  group.by = "celltype_to_plot", 
  slot = "scale.data"
)
human_agg_TAL <- AverageExpression(
  TAL_human_KTX_subclustering, 
  assays = "RNA", 
  features = human_genes_to_scale, 
  group.by = "celltype_to_plot",  
  slot = "scale.data"
)

gene_mapping_filtered_TAL <- gene_mapping_TAL[gene_mapping_TAL$MGI.symbol %in% rownames(mouse_agg_TAL[["RNA"]]) & 
                                                gene_mapping_TAL$HGNC.symbol %in% rownames(human_agg_TAL[["RNA"]]), ]

mouse_exp_TAL <- as.matrix(mouse_agg_TAL[["RNA"]][gene_mapping_filtered_TAL$MGI.symbol, ])
human_exp_TAL <- as.matrix(human_agg_TAL[["RNA"]][gene_mapping_filtered_TAL$HGNC.symbol, ])

cor_matrix_TAL <- cor(mouse_exp_TAL, human_exp_TAL, method = "spearman")


breaks_TAL = seq(0, 0.4, by = 0.02)

pheatmap(
  cor_matrix_TAL, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = FALSE, 
  color = colorRampPalette(c("white", "firebrick3"))(length(breaks_TAL) - 1),  # Diverging palette with white at 0
  breaks = breaks_TAL,  # More granular control
  legend = TRUE, 
  fontsize = 12,  
  fontsize_row = 12,  
  fontsize_col = 12   
)
