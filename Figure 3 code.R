#Figure 3

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)

setwd("D:/R_Ordner/KTx Daten")
Ktx_PT_mouse_subclustering <- readRDS("D:/R_Ordner/Ktx Daten/Ktx_PT_mouse_subclustering.rds")
Ktx_TAL_mouse_subclustering <- readRDS("D:/R_Ordner/Ktx Daten/Ktx_TAL_mouse_subclustering.rds")

#Figure 3A UMAP
Idents(Ktx_PT_mouse_subclustering) <- "celltype_level_2"
DimPlot(Ktx_PT_mouse_subclustering, label = F, cols = c('PT_S1' = '#CD5B45', 
                                        'PT_S2' = '#EE6A50', 
                                        'PT_S3' = '#8B3E2F', 
                                        'PT_S3_medullary' = '#EE7600' ,
                                        'PT_prolif' = '#8B668B' ,
                                        'PT_Injury_m1' = '#00BFFF' ,
                                        'PT_Injury_m2' = '#FFB90F' ,
                                        'PT_Injury_m3' = '#36648B' ,
                                        'PT_Injury_m4' = '#EE2C2C')) + NoLegend() 

#Figure 3A Heatmap
desired_order <- c("PT S1", "PT S2",  "PT S3", "PT S3 medullary", "PT prolif", "PT Injury m1", "PT Injury m2", "PT Injury m3", "PT Injury m4")

gene_order <- c(
  "Cd44", "Krt20", "Havcr1", "Klf6", "Vcam1", "Vim", 
  "Cxcl10", "Cxcl1", "Nfkbia", "Cdkn1a", "Ifit2", 
  "Gbp2", "Cd274", "Parp12", "Rnf19b", "Nqo1", 
  "Atad2", "Rad51b", "Pola1", "Mki67", "Top2a", 
  "Cyp7b1", "Slc22a7", "Slc7a13", "Slc22a6", "Slc34a1", "Slc7a7", "Slc5a2"
)

avg_exp <- AverageExpression(
  Ktx_PT_mouse_subclustering, 
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

# Figure 3B
df_mouse <- Ktx_PT_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2) %>%
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

celltype_order <- c('PT_S1', 'PT_S2', 'PT_S3', 'PT_S3_medullary',
                    'PT_prolif', 'PT_Injury_m1', 'PT_Injury_m2', 
                    'PT_Injury_m3', 'PT_Injury_m4')

df_mouse$celltype_level_2 <- factor(df_mouse$celltype_level_2, levels = celltype_order)

df_mouse_long <- df_mouse %>%
  select(ID_to_plot, celltype_level_2, percent) %>%
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

ggplot(df_mouse_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13)
  )

ggplot(df_mouse_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),    
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),  
    legend.position = "none"          
  )



#Figure 3C UMAP
Idents(Ktx_TAL_mouse_subclustering) <- "celltype_level_2"
DimPlot(Ktx_TAL_mouse_subclustering, label = F, cols = c('cTAL1' = '#4169E1', 
                                                           'cTAL2' = '#27408B', 
                                                           'mTAL1' = '#66CD00',
                                                           'mTAL2' = '#8FBC8F',
                                                           'cTAL3' = '#A4D3EE',
                                                           'cTAL4' = '#607B8B',
                                                           'TAL_Injury_m1' = '#EEC900',
                                                           'TAL_Injury_m2' = '#CD2626',
                                                           'TAL_Injury_m3' = '#AB82FF', 
                                                           'TAL_Prolif' = '#CD661D', 
                                                           'MD' = '#EEAEEE')) + NoLegend() 

#Figure 3C Heatmap
desired_order <- c("cTAL1", "cTAL2",  "cTAL3", "cTAL4", "mTAL1", "mTAL2", "MD", "TAL Prolif", "TAL Injury m1", "TAL Injury m2", "TAL Injury m3")

gene_order <- c(
  "Cd44", "Itga2", "Nfkbia", "Nfkbiz", "Tpm1", "Spp1", "Itgb1", 
  "Met", "Gbp2", "Gbp7", "Cxcl10", "Parp14", "Stat1", "Stat2", 
  "Irf1", "Iigp1", "Pole", "Atad2", "Top2a", "Robo2", "Pappa2", 
  "Nos1", "Cldn10", "Ranbp3l", "Sned1", "Tmem207", "Shd", 
  "Sfrp1", "Dusp15", "Tenm4", "Cldn16", "Kcnj10", "Enox1", "Tmem52b"
)

avg_exp <- AverageExpression(
  Ktx_TAL_mouse_subclustering, 
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


#Figure 3D 
df_tal <- Ktx_TAL_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(ID_to_plot) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_ids_tal <- c("BALB/c to BALB/c 1", "BALB/c to BALB/c 2", "C57BL/6 to C57BL/6 1",
                     "BALB/c to C57BL/6 1", "BALB/c to C57BL/6 2", 
                     "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3")

df_tal$ID_to_plot <- factor(df_tal$ID_to_plot, levels = ordered_ids_tal)

celltype_order_tal <- c("cTAL1", "cTAL2", "cTAL3", "cTAL4", "mTAL1", "mTAL2", "MD", 
                        "TAL_Prolif", "TAL_Injury_m1", 
                        "TAL_Injury_m2", "TAL_Injury_m3")

df_tal$celltype_level_2 <- factor(df_tal$celltype_level_2, levels = celltype_order_tal)

df_tal_long <- df_tal %>%
  select(ID_to_plot, celltype_level_2, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

ID_colors_tal <- c("BALB/c to BALB/c 1" = "grey70", 
                   "BALB/c to BALB/c 2" = "grey50", 
                   "C57BL/6 to C57BL/6 1" = "grey30",
                   "BALB/c to C57BL/6 1" = "darkolivegreen2",
                   "BALB/c to C57BL/6 2" = "darkolivegreen3",
                   "C57BL/6 to BALB/c 1" = "lightblue3",
                   "C57BL/6 to BALB/c 2" = "lightblue4",
                   "C57BL/6 to BALB/c 3" = "lightblue")

ggplot(df_tal_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_tal) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13)  
  )

ggplot(df_tal_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_tal) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),    
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),   
    legend.position = "none"          
  )

