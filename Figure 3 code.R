#Figure 3

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)


Ktx_PT_mouse_subclustering <- readRDS(".../Ktx_PT_mouse_subclustering.rds")
Ktx_TAL_mouse_subclustering <- readRDS(".../Ktx_TAL_mouse_subclustering.rds")

#Figure 3A UMAP

png("PT_mouse_umap.png", width = 400, height = 300, res = 200)  
DimPlot(
  Ktx_PT_mouse_subclustering,
  label = FALSE,
  cols = c(
    'PT_S1' = '#CD5B45', 
    'PT_S2' = '#EE6A50', 
    'PT_S3' = '#8B3E2F', 
    'PT_S3_medullary' = '#EE7600',
    'PT_prolif' = '#8B668B',
    'PT_Injury_m1' = '#00BFFF',
    'PT_Injury_m2' = '#FFB90F',
    'PT_Injury_m3' = '#36648B',
    'PT_Injury_m4' = '#EE2C2C'
  )
) + 
  theme_void() + 
  theme(legend.position = "none")
dev.off()

#for Figure 3A PAGA PT and TAL see the PT.PAGA.ipynb and TAL.PAGA.ipynb

#Figure 3C Heatmap
desired_order <- c("PT S1", "PT S2", "PT S3", "PT S3 med.", 
                   "PT prolif", "PT Injury m1", "PT Injury m2", 
                   "PT Injury m3", "PT Injury m4")

gene_order <- c(
  "Cd44", "Havcr1", "Vcam1",  "Anxa3",                                                      
  "Cxcl10", "Nfkb1",  "Atf3", 
  "Vim", "Itgb1", "Tnfaip3", "Tpm1", 
  "Birc3", "Rhob",     
  "Parp14", "Cd47", 
  "Nfkbia", "Icam1",           
  "Ptpre", "S100a10",          
  "Klf6", "Relb" ,  
  "Ifit2", "Ccl2",  
  "Klf7", "C3", "Pros1", "Sparc", "Cd74", "Stat1", "B2m", "Ppp2r2b", "Tspan8", "Pros1",  
  "Map3k1", "Map3k13", "Kcnh1", "Egfr" ,                                                     
  "Cryab", "Myo5b", "Cp", "Serpina10", "Jag1", "Wnt5a", "Foxc1", "Abcc1", "Gclc", "G6pdx", "Acsm3", "Stk39",       
  "Atad2", "Rad51b", "Pola1", "Mki67", "Top2a",                                            
  "Cyp7b1", "Nqo1", "Egf", "Slc7a12", "Slc7a13", "Slc22a6", "Lrp2", "Slc34a1", "Slc4a4", "Slc7a7", "Slc5a2" #S1-S3   #Lrp2
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

# highlight genes with "*", that are available in the spatial data                               
spat_genes <- scan("spat.genes.txt", what = character(), sep = "\n", quiet = TRUE)


row_labels <- ifelse(
  rownames(avg_exp_matrix_norm) %in% spat_genes,
  paste0(rownames(avg_exp_matrix_norm), "*"),
  rownames(avg_exp_matrix_norm)
)

pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = row_labels,   
  show_colnames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)

# Figure 3D
df_mouse <- Ktx_PT_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2_rev) %>%
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

df_mouse$celltype_level_2_rev <- factor(df_mouse$celltype_level_2_rev, levels = celltype_order)

df_mouse_long <- df_mouse %>%
  dplyr::select(ID_to_plot, celltype_level_2_rev, percent) %>%
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

ggplot(df_mouse_long, aes(x = celltype_level_2_rev, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13)
  )

ggplot(df_mouse_long, aes(x = celltype_level_2_rev, y = value, fill = ID_to_plot)) +
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

# Statistical testing with Benjamini Hochberg correction 
library(dplyr)
library(rstatix)
library(tidyr)

df_mouse <- Ktx_PT_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2_rev) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(ID_to_plot) %>%
  mutate(total = sum(count),
         percent = count / total) %>%
  ungroup()

group_map <- Ktx_PT_mouse_subclustering@meta.data %>%
  distinct(ID_to_plot, group)

df_mouse <- df_mouse %>%
  left_join(group_map, by = "ID_to_plot")

all_combinations <- expand.grid(
  ID_to_plot = unique(df_mouse$ID_to_plot),
  celltype_level_2_rev = unique(df_mouse$celltype_level_2_rev)
) %>%
  left_join(group_map, by = "ID_to_plot")

df_mouse_complete <- all_combinations %>%
  left_join(df_mouse, by = c("ID_to_plot", "celltype_level_2_rev", "group")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    total = ifelse(is.na(total), 0, total),
    percent = ifelse(is.na(percent), 0, percent)
  )

t_test_results <- df_mouse_complete %>%
  group_by(celltype_level_2_rev) %>%
  filter(length(unique(group)) == 2) %>%   
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%  
  select(celltype_level_2_rev, statistic, df, p, p.adj, p.adj.signif, group1, group2, n1, n2) %>%
  ungroup()

print(t_test_results)



#Figure 3F UMAP
Idents(Ktx_TAL_mouse_subclustering) <- "celltype_level_2"
DimPlot(Ktx_TAL_mouse_subclustering, label = F, cols = c('cTAL1' = '#4169E1', 
                                                           'cmTAL1' = '#27408B', 
                                                           'mTAL1' = '#66CD00',
                                                           'mTAL2' = '#8FBC8F',
                                                           'cmTAL2' = '#A4D3EE',
                                                           'cTAL2' = '#607B8B',
                                                           'TAL_Injury_m1' = '#EEC900',
                                                           'TAL_Injury_m2' = '#CD2626',
                                                           'TAL_Injury_m3' = '#AB82FF', 
                                                           'TAL_Prolif' = '#CD661D', 
                                                           'MD' = '#EEAEEE')) + NoLegend() 

#Figure 3C Heatmap
desired_order <- c("cTAL1", "cTAL2",  "cmTAL1", "cmTAL2", "mTAL1", "mTAL2", "MD", "TAL Prolif", "TAL Injury m1", "TAL Injury m2", "TAL Injury m3")

gene_order <- c(
 "Mmp7", "Klf6", "Spp1", "Atf3", "Pmepa1", "Fstl1", "Flna", "Pdgfb", "Fosl2", "Ets1", "Hbegf", "Lif", "Relb", "Tpm1", "Nfkbiz", "Nupr1",  "Serpina10",
  "F3", "Anxa2", "Met", "Mylk",   
  "Icam1", "Socs3", "Gbp2", "Cebpd", "Edn1", "Cxcl10", "Dcdc2a", "Sdc1", "Cp",                      
  "Parp9","Parp14", "Stat1", "Stat2", "Irf1", "Xaf1", "Epsti1", "Upp1" ,  "Tnfsf10",                             
  "Pola1", "Atad2", "Top2a", "Mki67",
  "Robo2", "Pappa2", "Nos1", "Bmp2", "Etv1", "Mgll", "Ccdc80", "Nox4", "Ltc4s", 
  "Cldn10", "Car2", "Ttn", "Slc2a1", "Prss23", "Clcnka", "Ppp1r1a", "Ppargc1a",  "Kcnj1", "Bmp5", "Selenbp1", "Kcnj10", "Tenm4", "Kynu",  "Kcnj16", "Cldn16", "Tmem52b" 
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

# highlight genes with "*", that are available in the spatial data    
spat_genes <- scan("spat.genes.txt", what = character(), sep = "\n", quiet = TRUE)

row_labels <- ifelse(rownames(avg_exp_matrix_norm) %in% spat_genes,
                     paste0("", rownames(avg_exp_matrix_norm), "*"),
                     rownames(avg_exp_matrix_norm))

pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = row_labels,
  show_colnames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100),
  border_color = NA,
  legend = FALSE
)

#Figure 3I 
library(dplyr)
library(tidyr)
library(ggplot2)

df_tal <- Ktx_TAL_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_to_plot) %>%
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

desired_order <- c("cTAL1", "cTAL2", "cmTAL1", "cmTAL2",
                   "mTAL1", "mTAL2", "MD",
                   "TAL Prolif", 
                   "TAL Injury m1", "TAL Injury m2", "TAL Injury m3")

df_tal$celltype_to_plot <- factor(df_tal$celltype_to_plot, levels = desired_order)

df_tal_long <- df_tal %>%
  select(ID_to_plot, celltype_to_plot, percent) %>%
  pivot_longer(cols = percent, names_to = "name", values_to = "value")

ID_colors_tal <- c("BALB/c to BALB/c 1" = "grey70", 
                   "BALB/c to BALB/c 2" = "grey50", 
                   "C57BL/6 to C57BL/6 1" = "grey30",
                   "BALB/c to C57BL/6 1" = "darkolivegreen2",
                   "BALB/c to C57BL/6 2" = "darkolivegreen3",
                   "C57BL/6 to BALB/c 1" = "lightblue3",
                   "C57BL/6 to BALB/c 2" = "lightblue4",
                   "C57BL/6 to BALB/c 3" = "lightblue")

ggplot(df_tal_long, aes(x = celltype_to_plot, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_tal) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 13))

ggplot(df_tal_long, aes(x = celltype_to_plot, y = value, fill = ID_to_plot)) +
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

# Statistical testing with Benjamini Hochberg correction                                
group_map <- Ktx_TAL_mouse_subclustering@meta.data %>%
  distinct(ID_to_plot, group)
df_tal <- df_tal %>%
  left_join(group_map, by = "ID_to_plot")

all_combinations <- expand.grid(
  ID_to_plot = unique(df_tal$ID_to_plot),
  celltype_to_plot = desired_order,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  left_join(group_map, by = "ID_to_plot")

df_tal_complete <- all_combinations %>%
  left_join(df_tal, by = c("ID_to_plot", "celltype_to_plot", "group")) %>%
  mutate(
    count   = ifelse(is.na(count), 0, count),
    total   = ifelse(is.na(total), 0, total),
    percent = ifelse(is.na(percent), 0, percent)
  )

ttest_results_tal <- df_tal_complete %>%
  group_by(celltype_to_plot) %>%
  filter(length(unique(group)) == 2) %>%             
  t_test(percent ~ group, var.equal = FALSE) %>%      
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(celltype_to_plot, statistic, df, p, p.adj, p.adj.signif,
         group1, group2, n1, n2) %>%
  ungroup()

print(ttest_results_tal, n = Inf)



