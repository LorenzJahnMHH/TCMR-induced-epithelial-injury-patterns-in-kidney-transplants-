# Figure 4
library(Seurat)
library(dplyr)
library(ggplot2)        
library(pheatmap)
library(dplyr)
library(rstatix)

Ktx_Leukos_mouse_subclustering <- readRDS(".../Ktx_Leukos_mouse_subclustering.rds")
# load mouse object for relative abundace of Leukocyte of all cells
Ktx_data_mouse <- readRDS(".../Ktx_data_mouse.rds") 

# Figure 4A UMAP
DimPlot(Ktx_Leukos_mouse_subclustering, group.by = 'celltype_level_3_renamed', label = F, 
        cols = c('Cd4 Th Cell' = '#666666', 
                 'Treg' = '#B3DFD6',
                 'Cd8 T Cell' = '#000000',
                 'Cycling Cd8 T Cell' = '#E5E5E5',
                 'Resident macrophage' = '#F4D166',
                 'Infiltrating macrophage' = '#B3C25B',
                 'Monocyte'= '#CAFF70',
                 'Act. macrophage'= '#5CA455',
                 'Cycling macrophage' = '#146C36',  
                 'NK' = '#68228B',
                 'B cell' = '#FFD39B',
                 'DC' = '#EEEE00'))  + NoLegend()


# Figure 4B barplot
leuko_celltypes <- unique(Ktx_Leukos_mouse_subclustering$celltype_level_3_renamed)

total_cells_by_sample <- Ktx_data_mouse@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(total_cells = n(), .groups = 'drop')

df_leukos <- Ktx_data_mouse@meta.data %>%
  filter(celltype_level_3_renamed %in% leuko_celltypes) %>%
  group_by(ID_to_plot, celltype_level_3_renamed) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(total_cells_by_sample, by = "ID_to_plot") %>%
  mutate(percent = count / total_cells) %>%
  ungroup()

all_combinations <- expand.grid(
  ID_to_plot = unique(df_leukos$ID_to_plot),
  celltype_level_3_renamed = leuko_celltypes
)

df_leukos <- all_combinations %>%
  left_join(df_leukos, by = c("ID_to_plot", "celltype_level_3_renamed")) %>%
  mutate(percent = ifelse(is.na(percent), 0, percent))

ordered_ids_leukos <- c(
  "BALB/c to BALB/c 1", "BALB/c to BALB/c 2", "C57BL/6 to C57BL/6 1",
  "BALB/c to C57BL/6 1", "BALB/c to C57BL/6 2", 
  "C57BL/6 to BALB/c 1", "C57BL/6 to BALB/c 2", "C57BL/6 to BALB/c 3"
)

df_leukos$ID_to_plot <- factor(df_leukos$ID_to_plot, levels = ordered_ids_leukos)

# 6. Farben definieren (wie gehabt)
ID_colors_leukos <- c(
  "BALB/c to BALB/c 1" = "grey70", 
  "BALB/c to BALB/c 2" = "grey50", 
  "C57BL/6 to C57BL/6 1" = "grey30",
  "BALB/c to C57BL/6 1" = "darkolivegreen2",
  "BALB/c to C57BL/6 2" = "darkolivegreen3",
  "C57BL/6 to BALB/c 1" = "lightblue3",
  "C57BL/6 to BALB/c 2" = "lightblue4",
  "C57BL/6 to BALB/c 3" = "lightblue"
)

ggplot(df_leukos, aes(x = celltype_level_3_renamed, y = percent, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13)
  ) +
  labs(x = "Cell Type Level 3 (renamed)", y = "Percentage", fill = "Mouse ID")

ggplot(df_leukos, aes(x = celltype_level_3_renamed, y = percent, fill = ID_to_plot)) +
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

ordered_celltypes <- c(
  "B cell",
  "DC",
  "Treg",
  "Cd4 Th Cell",
  "Cd8 T Cell",
  "Cycling Cd8 T Cell",
  "NK",
  "Monocyte",
  "Infiltrating macrophage",
  "Resident macrophage",
  "Act. macrophage",
  "Cycling macrophage"
)

df_leukos$celltype_level_3_renamed <- factor(
  df_leukos$celltype_level_3_renamed,
  levels = ordered_celltypes
)

ggplot(df_leukos, aes(x = celltype_level_3_renamed, y = percent, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_leukos) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13)
  ) +
  labs(x = "Cell Type Level 3 (renamed)", y = "Percentage", fill = "Mouse ID")


df_leukos$celltype_level_3_renamed <- factor(
  df_leukos$celltype_level_3_renamed,
  levels = ordered_celltypes
)

ggplot(df_leukos, aes(x = celltype_level_3_renamed, y = percent, fill = ID_to_plot)) +
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

# Statistical testing with Benjamini Hochberg correction 
total_cells_by_sample <- Ktx_data_mouse@meta.data %>%
  group_by(ID_to_plot) %>%
  summarise(total_cells = n(), .groups = "drop")

leuko_celltypes <- unique(Ktx_Leukos_mouse_subclustering$celltype_level_3_renamed)

df_leuko <- Ktx_data_mouse@meta.data %>%
  filter(celltype_level_3_renamed %in% leuko_celltypes) %>%
  group_by(ID_to_plot, celltype_level_3_renamed) %>%
  summarise(count = n(), .groups = "drop")

group_map <- Ktx_data_mouse@meta.data %>%
  distinct(ID_to_plot, group)

df_leuko <- df_leuko %>%
  left_join(total_cells_by_sample, by = "ID_to_plot") %>%
  left_join(group_map, by = "ID_to_plot") %>%
  mutate(percent = count / total_cells)  # Anteil an ALLEN Zellen

all_combinations <- expand.grid(
  ID_to_plot = unique(df_leuko$ID_to_plot),
  celltype_level_3_renamed = leuko_celltypes
) %>%
  left_join(group_map, by = "ID_to_plot") %>%
  left_join(total_cells_by_sample, by = "ID_to_plot")

df_leuko_complete <- all_combinations %>%
  left_join(df_leuko, by = c("ID_to_plot", "celltype_level_3_renamed", "group", "total_cells")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    percent = ifelse(is.na(percent), 0, percent)
  )

df_leuko_complete$group <- factor(df_leuko_complete$group, levels = c("syn", "alo"))

t_test_results_leuko <- df_leuko_complete %>%
  group_by(celltype_level_3_renamed) %>%
  filter(length(unique(group)) == 2) %>%
  t_test(percent ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(celltype_level_3_renamed, statistic, df, p, p.adj, p.adj.signif, group1, group2, n1, n2) %>%
  ungroup()

mean_diff <- df_leuko_complete %>%
  group_by(celltype_level_3_renamed, group) %>%
  summarise(mean_percent = mean(percent)*100, .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_percent) %>%
  mutate(diff_alo_minus_syn = alo - syn)

t_test_results_leuko <- t_test_results_leuko %>%
  left_join(mean_diff, by = "celltype_level_3_renamed")
print(t_test_results_leuko)

# Figure 4C Heatmap
Macrophage_markers_final <- c(
  "Ly6c2","Cebpb","Vcan","Chil3","S100a8","S100a9", 
  "Csf1r","Adgre1", # general
  "Apoe","Cd68","C1qa","C1qb","C1qc",
  "Spp1","Ccl12","Cxcl12", "Fabp4","Fabp5","Col1a1","Sparc","Trem2","Mrc1","Cx3cr1"
)

Final_markers <- c(
  "Cd3e","Cd3g","Cd4","Foxp3","Tnfrsf4",
  "Il7r","Tcf7","Cd8a","Cd8b1","Mki67","Top2a","Stmn1",
  "Cd79b","Cd19","Cd79a","Ms4a1",
  "Gzma","Klre1","Ncr1",
  "Flt3", "Xcr1","Clec9a","Plbd1",
  Macrophage_markers_final
)

avg_exp <- AverageExpression(
  Ktx_Leukos_mouse_subclustering, 
  features = Final_markers, 
  group.by = "celltype_level_3_renamed", 
  assays = "RNA", 
  slot = "data"
)

avg_exp_matrix <- avg_exp$RNA[Final_markers, ]
avg_exp_matrix_norm <- t(apply(avg_exp_matrix, 1, function(x) x / max(x, na.rm = TRUE)))


heatmap_plot <- pheatmap(
  avg_exp_matrix_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  breaks = seq(0, 1, length.out = 101),
  color = colorRampPalette(c("black", "yellow"))(100), # deine Farbskala
  border_color = NA,
  legend = FALSE
)



# Figure 4B
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

# Figure 4D - see folder domains_cooccurence/

