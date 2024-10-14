setwd("D:/R_Ordner/Ktx Daten")

#Figure 4A UMAP
Ktx_Leuko_mouse_subclustering <- readRDS("D:/R_Ordner/KTx Daten/Ktx_Leukos_mouse_subclustering_harmony.rds")

DimPlot(Ktx_Leuko_mouse_subclustering, group.by = 'celltype_level_2', label = F, 
        cols = c('Cd4+ T cells' = 'grey40', 
                 'Cd8+ T cells' = 'black',
                 'Cd8+ T cells_2_prolif' = 'grey90',
                 'Macro_1' = 'darkolivegreen',
                 'Macro_2'= 'darkolivegreen1',
                 'Macro_3'= 'darkolivegreen3',
                'Macro_prolif' = 'chartreuse4',  
                 'NK' = 'darkorchid4',
                 'B_cells' = 'burlywood1',
                 'DC' = 'yellow2'))  + NoLegend()

#Figure 4A Heatmap
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

celltype_order <- c("B_cells", "Cd4+ T cells", "Cd8+ T cells", "Cd8+ T cells_2_prolif", 
                    "DC", "Macro_1", "Macro_2", "Macro_3", "Macro_prolif", "NK")

gene_order <- c(
  "Klre1", "Gzma", "Ncr1", "Klrd1",        
  "C1qa", "C1qb", "C1qc",  "Col4a2", "Fstl1", "Col4a1", "Cd68", "Lyz2", "Fcer1g", "Bst1", "Ms4a7", "Apoe", "Ctss", "Adgre1", "Csf1r", "Creb5", "Cd300a",    # Macrophage markers
  "Naaa", "Plbd1", "Cbfa2t3", "Sept3",  
  "Top2a", "Mki67", "Kif11", "Diaph3",  
   "Cd8b1", "Cd8a", "Cd8b1", "Ms4a4b",        
  "Cd3e", "Cd4", "Tnfrsf4", "Rgs1", "Ltb",    
  "Ms4a1", "Cd79a", "Cd19", "Cd79b"      
)

cpm_matrix <- matrix(0, nrow = length(gene_order), ncol = length(celltype_order))
rownames(cpm_matrix) <- gene_order
colnames(cpm_matrix) <- celltype_order

for (ct in celltype_order) {
  cells.tmp <- rownames(Ktx_Leuko_mouse_subclustering@meta.data[
    Ktx_Leuko_mouse_subclustering$celltype_level_2 == ct, 
  ])
  
  if (length(cells.tmp) == 0) {
    next  
  }
  
  gene_counts <- Ktx_Leuko_mouse_subclustering@assays$RNA@counts[gene_order, cells.tmp, drop = FALSE]
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
  main = "Gene Expression in Leukocyte Subclusters",
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



#Figure 4C
total_cells_by_sample <- Ktx_data@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_leukos <- Ktx_Leuko_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample, by = "ID_to_plot") %>% 
  ungroup() %>%
  mutate(percent = count / total_cells)  

ordered_ids_leukos <- c("BALB/c to BALB/c 1", "BALB/c to BALB/c 2", "C57BL/6 to C57BL/6 1",
                        "BALB/c to C57BL/6 1", "BALB/c to C57BL/6 2", 
                        "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3")

df_leukos$ID_to_plot <- factor(df_leukos$ID_to_plot, levels = ordered_ids_leukos)

df_leukos_long <- df_leukos %>%
  select(ID_to_plot, celltype_level_2, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

ID_colors_leukos <- c("BALB/c to BALB/c 1" = "grey70", 
                      "BALB/c to BALB/c 2" = "grey50", 
                      "C57BL/6 to C57BL/6 1" = "grey30",
                      "BALB/c to C57BL/6 1" = "darkolivegreen2",
                      "BALB/c to C57BL/6 2" = "darkolivegreen3",
                      "C57BL/6 to BALB/c 1" = "lightblue3",
                      "C57BL/6 to BALB/c 2" = "lightblue4",
                      "C57BL/6 to BALB/c 3" = "lightblue")

ggplot(df_leukos_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

ggplot(df_leukos_long, aes(x = celltype_level_2, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Entferne die x-Achsenbeschriftungen
    axis.text.y = element_blank(),    # Entferne die y-Achsenbeschriftungen
    axis.title.x = element_blank(),   # Entferne den Titel der x-Achse
    axis.title.y = element_blank(),   # Entferne den Titel der y-Achse
    legend.position = "none",         # Entferne die Legende
    panel.grid.minor = element_blank() # Entferne Neben-Rasterlinien
  )


### with placeholders for syngeneic samples without cells

library(dplyr)
library(tidyr)
library(ggplot2)

total_cells_by_sample <- Ktx_data@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_leukos <- Ktx_Leuko_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample, by = "ID_to_plot") %>%
  mutate(percent = count / total_cells) %>%  # Calculate the relative abundance of cells
  ungroup()

all_combinations <- expand.grid(
  ID_to_plot = unique(df_leukos$ID_to_plot),
  celltype_level_2 = unique(df_leukos$celltype_level_2)
)

df_leukos <- all_combinations %>%
  left_join(df_leukos, by = c("ID_to_plot", "celltype_level_2")) %>%
  mutate(percent = ifelse(is.na(percent), 0, percent))  # Replace NA values with 0

ordered_ids_leukos <- c("BALB/c to BALB/c 1", "BALB/c to BALB/c 2", "C57BL/6 to C57BL/6 1",
                        "BALB/c to C57BL/6 1", "BALB/c to C57BL/6 2", 
                        "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3")

df_leukos$ID_to_plot <- factor(df_leukos$ID_to_plot, levels = ordered_ids_leukos)

ID_colors_leukos <- c("BALB/c to BALB/c 1" = "grey70", 
                      "BALB/c to BALB/c 2" = "grey50", 
                      "C57BL/6 to C57BL/6 1" = "grey30",
                      "BALB/c to C57BL/6 1" = "darkolivegreen2",
                      "BALB/c to C57BL/6 2" = "darkolivegreen3",
                      "C57BL/6 to BALB/c 1" = "lightblue3",
                      "C57BL/6 to BALB/c 2" = "lightblue4",
                      "C57BL/6 to BALB/c 3" = "lightblue")

ggplot(df_leukos, aes(x = celltype_level_2, y = percent, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  ) +
  labs(x = "Cell Type Level 2", y = "Percentage", fill = "Mouse ID")

ggplot(df_leukos, aes(x = celltype_level_2, y = percent, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),    
            axis.text.y = element_blank(),   
            axis.title.x = element_blank(),  
            axis.title.y = element_blank(),  
            legend.position = "none",         
            panel.grid.minor = element_blank())

#t test and benjamini hochberg 
library(dplyr)
library(rstatix)

total_cells_by_sample <- Ktx_data@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_leukos <- Ktx_Leuko_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_level_2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample, by = "ID_to_plot") %>%
  mutate(percent = count / total_cells) %>%
  ungroup() %>%
  mutate(group = Ktx_Leuko_mouse_subclustering$group[match(ID_to_plot, Ktx_Leuko_mouse_subclustering$ID_to_plot)])

unique_samples <- df_leukos %>%
  distinct(ID_to_plot, group)

all_combinations <- expand.grid(
  ID_to_plot = unique_samples$ID_to_plot,
  celltype_level_2 = unique(df_leukos$celltype_level_2)
) %>%
  left_join(unique_samples, by = "ID_to_plot")

df_leukos_complete <- all_combinations %>%
  left_join(df_leukos, by = c("ID_to_plot", "celltype_level_2", "group")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    total_cells = ifelse(is.na(total_cells), 0, total_cells),
    percent = ifelse(is.na(percent), 0, percent)
  )

t_test_results_leukos <- df_leukos_complete %>%
  group_by(celltype_level_2) %>%
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  select(celltype_level_2, statistic, df, p, p.adj, group1, group2, n1, n2) %>%
  ungroup()

t_test_results_leukos
