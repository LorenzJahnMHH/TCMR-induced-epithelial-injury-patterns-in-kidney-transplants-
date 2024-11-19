#Figure 4
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

setwd("D:/R_Ordner/Ktx Daten")
Ktx_Leuko_mouse_subclustering <- readRDS("D:/R_Ordner/KTx Daten/Ktx_Leukos_mouse_subclustering_harmony.rds")

#Figure 4A UMAP

DimPlot(Ktx_Leuko_mouse_subclustering, group.by = 'celltype_level_2', label = F, 
        cols = c('Cd4+ T cells' = '#666666', 
                 'Cd8+ T cells' = '#000000',
                 'Cd8+ T cells_2_prolif' = '#E5E5E5',
                 'Macro_1' = '#556B2F',
                 'Macro_2'= '#CAFF70',
                 'Macro_3'= '#A2CD5A',
                 'Macro_prolif' = '#458B00',  
                 'NK' = '#68228B',
                 'B_cells' = '#FFD39B',
                 'DC' = '#EEEE00'))  + NoLegend()

#Figure 4A Heatmap
celltype_order <- c("B_cells", "Cd4+ T cells", "Cd8+ T cells", "Cd8+ T cells_2_prolif", 
                    "DC", "Macro_1", "Macro_2", "Macro_3", "Macro_prolif", "NK")

celltype_order <- gsub("_", "-", celltype_order)

gene_order <- c(
  "Klre1", "Gzma", "Ncr1", "Klrd1",        
  "C1qa", "C1qb", "C1qc",  "Col4a2", "Fstl1", "Col4a1", "Cd68", "Lyz2", "Fcer1g", "Bst1", "Ms4a7", "Apoe", "Ctss", "Adgre1", "Csf1r", "Creb5", "Cd300a",   
  "Naaa", "Plbd1", "Cbfa2t3", "Sept3",  
  "Top2a", "Mki67", "Kif11", "Diaph3",  
  "Cd8b1", "Cd8a", "Ms4a4b",        
  "Cd3e", "Cd4", "Tnfrsf4", "Rgs1", "Ltb",    
  "Ms4a1", "Cd79a", "Cd19", "Cd79b"      
)

avg_exp <- AverageExpression(
  Ktx_Leuko_mouse_subclustering, 
  features = gene_order, 
  group.by = "celltype_level_2", 
  assays = "RNA", 
  slot = "data" 
)

avg_exp_matrix <- avg_exp$RNA[gene_order, celltype_order]
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



#Figure 4B
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
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),    
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),   
    legend.position = "none",         
    panel.grid.minor = element_blank() 
  )


# with placeholders for syngeneic samples without cells
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
  mutate(percent = count / total_cells) %>% 
  ungroup()

all_combinations <- expand.grid(
  ID_to_plot = unique(df_leukos$ID_to_plot),
  celltype_level_2 = unique(df_leukos$celltype_level_2)
)

df_leukos <- all_combinations %>%
  left_join(df_leukos, by = c("ID_to_plot", "celltype_level_2")) %>%
  mutate(percent = ifelse(is.na(percent), 0, percent))  

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

