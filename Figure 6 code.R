# Figure 6

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(biomaRt)

# Figure 6A UMAP
Ktx_data_human <- readRDS(".../Ktx_data_human.rds")

Idents(Ktx_data_human) <- "celltype_level_1"
DimPlot(Ktx_data_human, label = F, cols = c('PT' = 'coral1', 'tL' = 'brown4', 'TAL' = 'darkblue', 'PEC' = 'antiquewhite', 
                                                   'Leuko' = 'mediumpurple', 'Podo' = 'darkorange1', 
                                                   'IntC' = 'plum', 'EC' = 'lightskyblue1', 'CD-PC' = 'olivedrab3', 
                                                   'CNT' = 'olivedrab', 'DCT' = 'darkkhaki', 'Uro' = 'goldenrod1', 
                                                   'CD-IC' = 'ivory3', 'Myo' = 'ivory4'))    + NoLegend()

# Figure 6A Heatmap
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


#Figure 6B
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



# Figure 6C UMAP
Ktx_PT_human_subclustering <- readRDS(".../Ktx_PT_human_subclustering.rds")
Ktx_TAL_human_subclustering <- readRDS(".../Ktx_TAL_human_subclustering.rds")

DimPlot(Ktx_PT_human_subclustering, group.by = 'celltype_level_2', label = F, 
        cols = c('PT_S1' = 'coral3',
                 'PT_S2' = 'firebrick',  
                 'PT_S3' = 'coral4', 
                 'PT_S3_medullary' = 'darkorange2',
                 'PT_S1/2' = 'grey25', 
                 'PT_Injury_h1' = 'darkgoldenrod1', 
                 'PT_Injury_h2' = 'firebrick2'))  + NoLegend()

#Figure 5C Heatmap
desired_order <- c("PT S1", "PT S1/2", "PT S2",  "PT S3", "PT S3 med.", "PT Injury h1", "PT Injury h2")

gene_order <- rev(c("SLC5A2", "SLC7A7", "SLC4A4", "SLC34A1", "SLC22A6", "LRP2", "SLC7A13", "EGF",
                    "ACAA2", "DDT", "GHITM", "ETFB", "PHYH", "SLC1A5", "GABARAPL1", "PNP", "PGK1", "TXN", "FTL", "PRNP",
                    "ACSM3", "STK39", "MAP3K1", "KCNH1", "C3", "SPARC", "KLF7", "CCL2", "IFIT2", "RELB", "KLF6", "S100A10", "PTPRE", "ICAM1", "NFKBIA", "CD47",
                    "PARP14", "RHOB", "BIRC3", "TPM1", "TNFAIP3", "ITGB1", "VIM", "ATF3", "NFKB1", "CXCL10", "ANXA3", "VCAM1", "HAVCR1", "CD44"))

# manually annotate mouse markers to human heatmap
pt_m4 <- c("CCL2", "IFIT2", "RELB", "KLF6", "S100A10", "PTPRE", "ICAM1", "NFKBIA", "CD47",
           "PARP14", "RHOB", "BIRC3", "TPM1", "TNFAIP3", "ITGB1", "VIM", "ATF3", "NFKB1", "CXCL10", 
           "ANXA3", "VCAM1", "HAVCR1", "CD44")
pt_m3 <- c("C3", "SPARC", "KLF7")
pt_m2 <- c("MAP3K1", "KCNH1")
pt_m1 <- c("ACSM3", "STK39")

avg_exp <- AverageExpression(
  Ktx_PT_human_subclustering, 
  features = gene_order, 
  group.by = "celltype_to_plot", 
  assays = "RNA", 
  slot = "data"
)

avg_exp_matrix <- avg_exp$RNA[gene_order, desired_order]
avg_exp_matrix_norm <- t(apply(avg_exp_matrix, 1, function(x) x / max(x, na.rm = TRUE)))

row_labels_pt <- sapply(rownames(avg_exp_matrix_norm), function(g) {
  if (g %in% pt_m4) return(paste0(g, " (m4)"))
  if (g %in% pt_m3) return(paste0(g, " (m3)"))
  if (g %in% pt_m2) return(paste0(g, " (m2)"))
  if (g %in% pt_m1) return(paste0(g, " (m1)"))
  return(g)
})
pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_row = row_labels_pt,
  show_colnames = TRUE,
  show_rownames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)

                               pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none", 
  main = "Marker Genes in PT Human Subclusters (Avg. Expression)",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA
)

                        
#Figure 6D
df_human <- Ktx_PT_human_subclustering@meta.data %>%
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

# Statistical testing with Benjamini Hochberg correction 
  total_cells_by_sample_human <- Ktx_PT_human_subclustering @meta.data %>%
  group_by(sample) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_human <- Ktx_PT_human_subclustering @meta.data %>%
  group_by(sample, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample_human, by = "sample") %>%
  mutate(percent = count / total_cells) %>%
  ungroup() %>%
  mutate(group = Ktx_PT_human_subclustering $group[match(sample, Ktx_PT_human_subclustering $sample)])

unique_samples_human <- df_human %>%
  distinct(sample, group)

all_combinations_human <- expand.grid(
  sample = unique_samples_human$sample,
  celltype_level_2 = unique(df_human$celltype_level_2)
) %>%
  left_join(unique_samples_human, by = "sample")

df_human_complete <- all_combinations_human %>%
  left_join(df_human, by = c("sample", "celltype_level_2", "group")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    total_cells = ifelse(is.na(total_cells), 0, total_cells),
    percent = ifelse(is.na(percent), 0, percent)
  )

t_test_results_human <- df_human_complete %>%
  group_by(celltype_level_2) %>%
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  select(celltype_level_2, statistic, df, p, p.adj, group1, group2, n1, n2) %>%
  ungroup()

print(t_test_results_human)

#Figure 6E 
Ktx_PT_mouse_subclustering <- readRDS("D:/R_Ordner/KTx Daten/Ktx Objekte final/Maus/Ktx_PT_mouse_subclustering.rds")
Ktx_PT_human_subclustering <- readRDS("D:/R_Ordner/KTx Daten/Ktx Objekte final/Human/Ktx_PT_human_subclustering.rds")

mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                      host = "https://dec2021.archive.ensembl.org/")
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                      host = "https://dec2021.archive.ensembl.org/")

gene_map <- getLDS(
  attributes = c("mgi_symbol"),
  filters = "mgi_symbol",
  values = rownames(Ktx_PT_mouse_subclustering),
  mart = mouse_mart,
  attributesL = c("hgnc_symbol"),
  martL = human_mart,
  uniqueRows = TRUE
)
gene_map <- gene_map[!duplicated(gene_map[, 2]), ]
colnames(gene_map) <- c("mouse_gene", "human_gene")

mouse_expr <- GetAssayData(Ktx_PT_mouse_subclustering, slot = "data")
mouse_expr <- mouse_expr[gene_map$mouse_gene, ]
rownames(mouse_expr) <- gene_map$human_gene

human_expr <- GetAssayData(Ktx_PT_human_subclustering, slot = "data")
common_genes <- intersect(rownames(mouse_expr), rownames(human_expr))
mouse_expr <- mouse_expr[common_genes, ]
human_expr <- human_expr[common_genes, ]


colnames(human_expr) <- make.unique(colnames(human_expr))
colnames(Ktx_PT_human_subclustering) <- make.unique(colnames(Ktx_PT_human_subclustering))


ref_labels_human <- Ktx_PT_human_subclustering@meta.data[colnames(human_expr), "celltype_level_2"]


singleR_results_PT <- SingleR(test = as.matrix(mouse_expr),
                              ref = as.matrix(human_expr),
                              labels = ref_labels_human)

Ktx_PT_mouse_subclustering$SingleR.labels.human <- singleR_results_PT$pruned.labels

confusion_PT <- with(
  Ktx_PT_mouse_subclustering@meta.data,
  table(celltype_level_2_rev, SingleR.labels.human)
)

prop_PT_mouse <- prop.table(confusion_PT, margin = 1)

pretty_rows <- c(
  "PT_S1"           = "PT S1",
  "PT_S1/2"         = "PT S1/2",
  "PT_S2"           = "PT S2",
  "PT_S3"           = "PT S3",
  "PT_S3_medullary" = "PT S3 med.",
  "PT_prolif"       = "PT prolif",
  "PT_Injury_m1"    = "PT Injury m1",
  "PT_Injury_m2"    = "PT Injury m2",
  "PT_Injury_m3"    = "PT Injury m3",
  "PT_Injury_m4"    = "PT Injury m4"
)
rownames(prop_PT_mouse) <- pretty_rows[rownames(prop_PT_mouse)]

pretty_cols <- c(
  "PT_S1"           = "PT S1",
  "PT_S1/2"         = "PT S1/2",
  "PT_S2"           = "PT S2",
  "PT_S3"           = "PT S3",
  "PT_S3_medullary" = "PT S3 med.",
  "PT_prolif"       = "PT prolif",
  "PT_Injury_h1"    = "PT Injury h1",
  "PT_Injury_h2"    = "PT Injury h2"
)
colnames(prop_PT_mouse) <- pretty_cols[colnames(prop_PT_mouse)]

pheatmap(
  prop_PT_mouse,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color        = colorRampPalette(c("white", "firebrick3"))(100),
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col    = 45,
  main         = "Mouse PT clusters mapped to human PT reference"
)



#Figure 6F UMAP
Idents(Ktx_TAL_mouse_subclustering ) <- "celltype_level_2"
DimPlot(Ktx_TAL_mouse_subclustering , group.by = 'celltype_level_2', label = F, 
        cols = c('cTAL1' = 'royalblue', 
                 'cTAL2' = 'royalblue4',
                 'mTAL1' = 'chartreuse3', 
                 'mTAL2' = 'darkseagreen',
                 'mTAL3' = 'darkgreen', 
                 'cTAL3' = 'lightskyblue2', 
                 'cTAL4' = 'lightskyblue4', 
                 'TAL_Injury_h1' = 'mediumpurple1',
                 'TAL_Injury_h2' = 'firebrick3')) + NoLegend()  

#Figure 6F Heatmap
desired_order_tal <- c("cTAL1", "cTAL2", "cTAL3", "cTAL4", "mTAL1", "mTAL2", "mTAL3", "TAL Injury h1", "TAL Injury h2")

gene_order_tal <- rev(c("TMEM52B", "CLDN16",  "KCNJ10", "TENM4", "HCRTR2",  "CABP1", "CALCR", "NR2E3", "PAG1" , "LYPD6", "KCNJ16", "KCNJ1", "MEST", "CLCNKA", "RDH8", "FAM43A", "NEFL", "SLC22A2","CNTD2",  "GAMT", "PPP1R1A",
                        "IGFBP6", "CCND2", "ENO2", "KLF2", "KLF4", "PDLIM1", "PRDX1", "PRDX4", "PRDX6", "RHOB", "JUN",
                        "PARP14", "STAT1", "IRF1", "SOCS3", "EDN1", "CEBPD", "F3", "NFKBIZ", "TPM1", "RELB", "LIF", "HBEGF", "ETS1", "FOSL2", "PDGFB",
                        "FLNA", "FSTL1", "PMEPA1", "ATF3", "SPP1", "MMP7"))

tal_m3 <- c("F3", "NFKBIZ", "TPM1", "RELB", "LIF", "HBEGF", "ETS1", "FOSL2", "PDGFB",
            "FLNA", "FSTL1", "PMEPA1", "ATF3", "SPP1", "MMP7")
tal_m2 <- c("SOCS3", "EDN1", "CEBPD")
tal_m1 <- c("PARP14", "STAT1", "IRF1")

avg_exp_tal <- AverageExpression(
  Ktx_TAL_human_subclustering, 
  features = gene_order_tal, 
  group.by = "celltype_to_plot", 
  assays = "RNA", 
  slot = "data"
)

avg_exp_matrix_tal <- avg_exp_tal$RNA[gene_order_tal, desired_order_tal]
avg_exp_matrix_tal_norm <- t(apply(avg_exp_matrix_tal, 1, function(x) x / max(x, na.rm = TRUE)))

row_labels_tal <- sapply(rownames(avg_exp_matrix_tal_norm), function(g) {
  if (g %in% tal_m3) return(paste0(g, " (m3)"))
  if (g %in% tal_m2) return(paste0(g, " (m2)"))
  if (g %in% tal_m1) return(paste0(g, " (m1)"))
  return(g)
})

pheatmap(
  avg_exp_matrix_tal_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_row = row_labels_tal,
  show_colnames = TRUE,
  show_rownames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)


#Figure 6G
library(Seurat)
library(dplyr)
library(ggplot2)

desired_order <- c('cTAL1', 'cTAL2', 'cTAL3', 'cTAL4', 
                   'mTAL1', 'mTAL2', 'mTAL3', 
                   'TAL_Injury_h1', 'TAL_Injury_h2')

df_tal <- Ktx_TAL_mouse_subclustering@meta.data %>%
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



#Figure 6H 
Ktx_TAL_mouse_subclustering <- readRDS(".../Ktx_TAL_mouse_subclustering.rds")
Ktx_TAL_human_subclustering <- readRDS(".../Human/Ktx_TAL_human_subclustering.rds")

gene_map <- getLDS(
  attributes = c("mgi_symbol"),
  filters = "mgi_symbol",
  values = rownames(Ktx_TAL_mouse_subclustering),
  mart = mouse_mart,
  attributesL = c("hgnc_symbol"),
  martL = human_mart,
  uniqueRows = TRUE
)
gene_map <- gene_map[!duplicated(gene_map[, 2]), ]
colnames(gene_map) <- c("mouse_gene", "human_gene")

mouse_expr <- GetAssayData(Ktx_TAL_mouse_subclustering, slot = "data")
mouse_expr <- mouse_expr[gene_map$mouse_gene, ]
rownames(mouse_expr) <- gene_map$human_gene

human_expr <- GetAssayData(Ktx_TAL_human_subclustering, slot = "data")
common_genes <- intersect(rownames(mouse_expr), rownames(human_expr))
mouse_expr <- mouse_expr[common_genes, ]
human_expr <- human_expr[common_genes, ]

colnames(human_expr) <- make.unique(colnames(human_expr))
colnames(Ktx_TAL_human_subclustering) <- make.unique(colnames(Ktx_TAL_human_subclustering))

ref_labels_human <- Ktx_TAL_human_subclustering@meta.data[colnames(human_expr), "celltype_level_2"]

singleR_results_TAL <- SingleR(test = as.matrix(mouse_expr),
                               ref = as.matrix(human_expr),
                               labels = ref_labels_human)

Ktx_TAL_mouse_subclustering$SingleR.labels.human_pruned <- singleR_results_TAL$pruned.labels

confusion_TAL <- with(
  Ktx_TAL_mouse_subclustering@meta.data,
  table(celltype_level_2, SingleR.labels.human_pruned)
)

row_order_tal_mouse <- c(
  "cTAL1",
  "cTAL2",
  "cmTAL1",
  "cmTAL2",
  "mTAL1",
  "mTAL2",
  "MD",
  "TAL_Prolif",
  "TAL_Injury_m1",
  "TAL_Injury_m2",
  "TAL_Injury_m3"
)

row_keep <- row_order_tal_mouse[row_order_tal_mouse %in% rownames(confusion_TAL)]
confusion_TAL <- confusion_TAL[row_keep, , drop = FALSE]

pretty_colnames_tal <- c(
  "cTAL1"         = "cTAL1",
  "cTAL2"         = "cTAL2",
  "cTAL3"         = "cTAL3",
  "cTAL4"         = "cTAL4",
  "mTAL1"         = "mTAL1",
  "mTAL2"         = "mTAL2",
  "mTAL3"         = "mTAL3",
  "TAL_Injury_h1" = "TAL Injury h1",
  "TAL_Injury_h2" = "TAL Injury h2"
)

colnames(confusion_TAL) <- pretty_colnames_tal[colnames(confusion_TAL)]

col_order_human_tal <- c(
  "cTAL1",
  "cTAL2",
  "cTAL3",
  "cTAL4",
  "mTAL1",
  "mTAL2",
  "mTAL3",
  "TAL Injury h1",
  "TAL Injury h2"
)

col_keep <- col_order_human_tal[col_order_human_tal %in% colnames(confusion_TAL)]
confusion_TAL <- confusion_TAL[, col_keep, drop = FALSE]

prop_TAL_mouse <- prop.table(confusion_TAL, margin = 1)
prop_TAL_mouse[is.na(prop_TAL_mouse)] <- 0

pheatmap(
  prop_TAL_mouse,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color        = colorRampPalette(c("white", "firebrick3"))(100),
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col    = 45,
  main         = "Mouse → Human Mapping (TAL)"
)


