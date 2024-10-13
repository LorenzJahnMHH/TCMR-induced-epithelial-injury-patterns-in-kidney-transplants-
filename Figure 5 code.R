# Figure 5

# Figure 5A UMAP
Ktx_human_int_harmony <- readRDS("Ktx_human_int_harmony.rds")

Idents(Ktx_human_int_harmony) <- "celltype_level_1"
DimPlot(Ktx_human_int_harmony, label = F, cols = c('PT' = 'coral1', 'tL' = 'brown4', 'TAL' = 'darkblue', 'PEC' = 'antiquewhite', 
                                                   'Leuko' = 'mediumpurple', 'Podo' = 'darkorange1', 
                                                   'IntC' = 'plum', 'EC' = 'lightskyblue1', 'CD-PC' = 'olivedrab3', 
                                                   'CNT' = 'olivedrab', 'DCT' = 'darkkhaki', 'Uro' = 'goldenrod1', 
                                                   'CD-IC' = 'ivory3', 'Myo' = 'ivory4'))    + NoLegend()

# Figure 5A Heatmap
Idents(Ktx_human_int_harmony) <- "celltype_level_1"

top_markers <- FindAllMarkers(
  object = Ktx_human_int_harmony, 
  only.pos = TRUE, 
  min.pct = 0.5,
  logfc.threshold = 0.5  
)

celltype_order = c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", 
                   "CD-IC",  "EC", "Leuko", "IntC", 
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

cpm_matrix <- matrix(0, nrow = length(top_genes_order), ncol = length(celltype_order))
rownames(cpm_matrix) <- top_genes_order
colnames(cpm_matrix) <- celltype_order

for (ct in celltype_order) {
  cells.tmp <- rownames(Ktx_human_int_harmony@meta.data[Ktx_human_int_harmony$celltype_level_1 == ct, ])
  
  if (length(cells.tmp) == 0) {
    next  
  }
  
  gene_counts <- Ktx_human_int_harmony@assays$RNA@counts[top_genes_order, cells.tmp, drop = FALSE]
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

#Figure 5B
df_human <- Ktx_human_int_harmony@meta.data %>%
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

DimPlot(PT_human_KTX_subclustering, group.by = 'celltype_updated', label = F, 
        cols = c('PT_S1' = 'coral3',
                 'PT_S2' = 'firebrick',  
                 'PT_S3' = 'coral4', 
                 'PT_S3_medullary' = 'darkorange2',
                 'PT_S1/2' = 'grey25', 
                 'PT_Injury_1' = 'firebrick2', 
                 'PT_Injury_2' = 'darkgoldenrod1'))  + NoLegend()

#Figure 5C Heatmap
desired_order <- c("PT S1", "PT S1/2", "PT S2",  "PT S3", "PT S3 med.", "PT Injury h1", "PT Injury h2")

gene_order <- rev(c("SLC5A2", "SLC7A7", "SLC22A8", "SLC34A1", "SLC7A13", "SLC22A24", "SLC5A1", "MT-ND1", "MT-ND2", "UQCRQ", "GPX4", "ACAA2", "ATP5F1E", "COX5B",
                    "IFITM3", "VIM", "VCAM1", "TMSB10", "HLA-A"))

cpm_matrix <- matrix(0, nrow = length(gene_order), ncol = length(desired_order))
rownames(cpm_matrix) <- gene_order
colnames(cpm_matrix) <- desired_order

for (ct in desired_order) {
  cells.tmp <- rownames(PT_human_KTX_subclustering@meta.data[PT_human_KTX_subclustering$celltype_to_plot == ct, ])
  
  if (length(cells.tmp) == 0) {
    next  
  }
  
  gene_counts <- PT_human_KTX_subclustering@assays$RNA@counts[gene_order, cells.tmp, drop = FALSE]
  gene_counts <- as.matrix(gene_counts)
  
  total_counts <- colSums(gene_counts)
  total_counts[total_counts == 0] <- 1
  
  cpm <- t(t(gene_counts) / total_counts * 1e6)
  cpm_matrix[, ct] <- rowMeans(cpm, na.rm = TRUE)
}

cpm_matrix_norm <- t(apply(cpm_matrix, 1, function(x) x / max(x, na.rm = TRUE)))

pheatmap(
  cpm_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  main = "Marker Genes in PT Human Subclusters",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA
)

pheatmap(
  cpm_matrix_norm,
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


#Figure 7D
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

df_human <- PT_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_samples <- c("C1", "C2", "C4", "R1", "R3", "R4")

df_human$sample <- factor(df_human$sample, levels = ordered_samples)

celltype_order <- c("PT_S1", "PT_S2", "PT_S1/2", "PT_S3", "PT_S3_medullary", "PT_Injury_1", "PT_Injury_2")

df_human$celltype_updated <- factor(df_human$celltype_updated, levels = celltype_order)

df_human_long <- df_human %>%
  select(sample, celltype_updated, percent) %>%
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

ggplot(df_human_long, aes(x = celltype_updated, y = value, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 13) 
  )

ggplot(df_human_long, aes(x = celltype_updated, y = value, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Remove x-axis labels
    axis.text.y = element_blank(),    # Remove y-axis labels
    axis.title.x = element_blank(),   # Remove x-axis title
    axis.title.y = element_blank(),   # Remove y-axis title
    legend.position = "none"          # Remove legend
  )

# t- test and Benjamini Hochberg
library(dplyr)
library(rstatix)
total_cells_by_sample_human <- PT_human_KTX_subclustering@meta.data %>%
  group_by(sample) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_human <- PT_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample_human, by = "sample") %>%
  mutate(percent = count / total_cells) %>%
  ungroup() %>%
  mutate(group = PT_human_KTX_subclustering$group[match(sample, PT_human_KTX_subclustering$sample)])

unique_samples_human <- df_human %>%
  distinct(sample, group)

all_combinations_human <- expand.grid(
  sample = unique_samples_human$sample,
  celltype_updated = unique(df_human$celltype_updated)
) %>%
  left_join(unique_samples_human, by = "sample")

df_human_complete <- all_combinations_human %>%
  left_join(df_human, by = c("sample", "celltype_updated", "group")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    total_cells = ifelse(is.na(total_cells), 0, total_cells),
    percent = ifelse(is.na(percent), 0, percent)
  )

t_test_results_human <- df_human_complete %>%
  group_by(celltype_updated) %>%
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  select(celltype_updated, statistic, df, p, p.adj, group1, group2, n1, n2) %>%
  ungroup()

print(t_test_results_human)

#Figure 7D Correlation heatmap
library(Seurat)
library(dplyr)
library(biomaRt)
library(pheatmap)

Idents(PT_human_KTX_subclustering) <- "celltype_to_plot"  
Idents(Ktx_PT_mouse_subclustering) <- "celltype_to_plot"  

PT_mouse_marker <- FindAllMarkers(
  Ktx_PT_mouse_subclustering, 
  only.pos = TRUE, 
  logfc.threshold = 1
)
PT_mouse_marker <- PT_mouse_marker %>%
  filter(p_val_adj < 0.05, pct.1 > 0.25)

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


pheatmap(
  cor_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = FALSE, 
  color = colorRampPalette(c("white", "firebrick3"))(100),
  legend = FALSE, 
  fontsize = 12,  
  fontsize_row = 12,  
  fontsize_col = 12   
)


#Figure 7E UMAP
Idents(TAL_human_KTX_subclustering ) <- "celltype_updated"
DimPlot(TAL_human_KTX_subclustering , group.by = 'celltype_updated', label = F, 
        cols = c('cTAL1' = 'royalblue', 
                 'cTAL2' = 'royalblue4',
                 'mTAL1' = 'chartreuse3', 
                 'mTAL2' = 'darkseagreen',
                 'mTAL3' = 'darkgreen', 
                 'cTAL3' = 'lightskyblue2', 
                 'cTAL4' = 'lightskyblue4', 
                 'TAL_Injury_1' = 'firebrick3',
                 'TAL_Injury_2' = 'mediumpurple1')) + NoLegend()  

#Figure 7E Heatmap
library(Seurat)
library(dplyr)
library(pheatmap)

desired_order_tal <- c("cTAL1", "cTAL2", "cTAL3", "cTAL4", "mTAL1", "mTAL2", "mTAL3", "TAL Injury h1", "TAL Injury h2")

gene_order_tal <- rev(c(

  "CLDN16", "TMEM52B", "KCNJ10",
  "PHACTR1", "CABP1", "TMEM72-AS1",
  "CALCR", "PAPPA2", "NOS1",
  "FGL2", "PDE3A", "PDE3A",  
  "RDH8", "PCBP3", "FMN1",
  "ATP5MC3", "COX7B", "SLC25A5",
  "MT-CYB", "MT-CO2", "MT-ND4",
  "KLF2", "CRIP1", "NOV", "IGFBP6", "SLPI",
  "SPP1", "SAMD4A", "ADAMTS1", "CD44", "MMP7"

))

cpm_matrix_tal <- matrix(0, nrow = length(gene_order_tal), ncol = length(desired_order_tal))
rownames(cpm_matrix_tal) <- gene_order_tal
colnames(cpm_matrix_tal) <- desired_order_tal

for (ct in desired_order_tal) {
  cells.tmp <- rownames(TAL_human_KTX_subclustering@meta.data[TAL_human_KTX_subclustering$celltype_to_plot == ct, ])
  
  if (length(cells.tmp) == 0) {
    next  
  }
  
  gene_counts <- TAL_human_KTX_subclustering@assays$RNA@counts[gene_order_tal, cells.tmp, drop = FALSE]
  gene_counts <- as.matrix(gene_counts)
  
  total_counts <- colSums(gene_counts)
  total_counts[total_counts == 0] <- 1
  
  cpm <- t(t(gene_counts) / total_counts * 1e6)
  cpm_matrix_tal[, ct] <- rowMeans(cpm, na.rm = TRUE)
}

cpm_matrix_tal_norm <- t(apply(cpm_matrix_tal, 1, function(x) x / max(x, na.rm = TRUE)))

pheatmap(
  cpm_matrix_tal_norm,
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


#Figure 7F
library(Seurat)
library(dplyr)
library(ggplot2)

desired_order <- c('cTAL1', 'cTAL2', 'cTAL3', 'cTAL4', 
                   'mTAL1', 'mTAL2', 'mTAL3', 
                   'TAL_Injury_h1', 'TAL_Injury_h2')

df_tal <- TAL_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

df_tal$celltype_updated <- factor(df_tal$celltype_updated, levels = desired_order)

sample_colors <- c("C1" = "grey70", 
                   "C2" = "grey50", 
                   "C4" = "grey30",
                   "R1" = "lightcoral",
                   "R3" = "indianred",
                   "R4" = "firebrick")

ggplot(df_tal, aes(x = celltype_updated, y = percent, fill = sample)) +
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

ggplot(df_tal, aes(x = celltype_updated, y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Entferne die x-Achsenbeschriftungen
    axis.text.y = element_blank(),    # Entferne die y-Achsenbeschriftungen
    axis.title.x = element_blank(),   # Entferne den Titel der x-Achse
    axis.title.y = element_blank(),   # Entferne den Titel der y-Achse
    legend.position = "none"          # Entferne die Legende
  )

library(Seurat)
library(dplyr)
library(rstatix)
library(ggplot2)

df_tal <- TAL_human_KTX_subclustering@meta.data %>%
  group_by(sample, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

desired_order <- c('cTAL1', 'cTAL2', 'cTAL3', 'cTAL4', 
                   'mTAL1', 'mTAL2', 'mTAL3', 
                   'TAL_Injury_h1', 'TAL_Injury_h2')

df_tal$celltype_updated <- factor(df_tal$celltype_updated, levels = desired_order)

df_tal <- df_tal %>%
  mutate(group = ifelse(sample %in% c("C1", "C2", "C4"), "stable.allo", "tcmr"))

t_test_results_tal <- df_tal %>%
  group_by(celltype_updated) %>%
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  select(celltype_updated, statistic, df, p, p.adj, group1, group2, n1, n2) %>%
  ungroup()

print(t_test_results_tal %>%
        mutate(p_value_before_correction = p, 
               p_value_after_correction = p.adj) %>%
        select(celltype_updated, n1, n2, p_value_before_correction, p_value_after_correction))


#Figure 7F Correlation heatmap
library(Seurat)
library(dplyr)
library(biomaRt)
library(pheatmap)

Idents(TAL_human_KTX_subclustering) <- "celltype_to_plot"  
Idents(Ktx_TAL_mouse_subclustering) <- "celltype_to_plot" 

TAL_marker <- FindAllMarkers(
  Ktx_TAL_mouse_subclustering, 
  only.pos = TRUE, 
  logfc.threshold = 1
)
TAL_marker <- TAL_marker %>%
  filter(p_val_adj < 0.05, pct.1 > 0.25)

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


pheatmap(
  cor_matrix_TAL, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = FALSE, 
  color = colorRampPalette(c("white", "firebrick3"))(100),
  legend = FALSE, 
  fontsize = 12, 
  fontsize_row = 12, 
  fontsize_col = 12  
)
