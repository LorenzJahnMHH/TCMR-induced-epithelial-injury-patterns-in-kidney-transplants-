# Figure 1

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)

Ktx_data_mouse_mouse <- readRDS(".../Ktx_data_mouse_mouse.rds") 


# Figure 1C UMAP
Idents(Ktx_data_mouse) <- "celltype_level_1"

DimPlot(Ktx_data_mouse, label = F, cols = c('PT' = '#FF7256', 'tL' = '#8B2323', 'TAL' = '#00008B', 'PEC' = '#FAEBD7', 
                                      'Leuko' = '#9370DB', 'Podo' = '#FF7F00', 
                                      'IntC' = '#DDA0DD', 'EC' = '#B0E2FF', 'CD-PC' = '#9ACD32', 
                                      'CNT' = '#6B8E23', 'DCT' = '#BDB76B', 'Prolif' = '#00C5CD', 'Uro' = '#FFC125', 
                                      'CD-IC-A' = '#CDCDC1', 'CD-IC-B' = '#8B8B83')) + NoLegend() 

#Figure 1C Heatmap
Idents(Ktx_data_mouse) <- "celltype_level_1"

top_markers <- FindAllMarkers(
  object = Ktx_data_mouse, 
  only.pos = TRUE, 
  min.pct = 0.25,
  logfc.threshold = 0.5  
)

celltype_order = c("Podo", "PT", "tL", "TAL", "DCT", "CNT", "CD-PC", 
                   "CD-IC-A", "CD-IC-B", "EC", "Leuko", "IntC", 
                   "Uro", "PEC", "Prolif")

top70 <- top_markers %>%
  group_by(cluster) %>%
  slice_head(n = 70) %>%  
  ungroup()


top_genes <- top70$gene

top_genes_order <- top70 %>% 
  mutate(cluster = factor(cluster, levels = rev(celltype_order))) %>%  
  arrange(cluster) %>% 
  pull(gene)

top_genes_order <- rev(top_genes_order)

avg_expression <- AverageExpression(
  Ktx_data_mouse,
  features = top_genes,
  group.by = "celltype_level_1",
  assays = "RNA", slot = "data"
)$RNA

avg_expression <- avg_expression[top_genes_order, celltype_order, drop = FALSE]

normalized_matrix <- t(apply(avg_expression, 1, function(x) x / max(x, na.rm = TRUE)))

heatmap_df <- as.data.frame(as.table(normalized_matrix))
colnames(heatmap_df) <- c("Gene", "CellType", "Expression")

heatmap_df$Gene <- factor(heatmap_df$Gene, levels = unique(top_genes_order))
heatmap_df$CellType <- factor(heatmap_df$CellType, levels = celltype_order)

p1 <- ggplot(heatmap_df, aes(x = CellType, y = Gene, fill = Expression)) + 
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

print(p1)

p2 <- ggplot(heatmap_df, aes(x = CellType, y = Gene, fill = Expression)) + 
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

print(p2)

# Figure 1E
df_mouse <- Ktx_data_mouse@meta.data %>%
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

celltype_order <- c("Podo", "PT", "tL", "TAL", "MD", "DCT", "CNT", 
                    "CD-PC", "CD-IC-A", "CD-IC-B", "EC", "Leuko", 
                    "IntC", "Uro", "PEC", "Prolif")

df_mouse$celltype_level_1 <- factor(df_mouse$celltype_level_1, levels = celltype_order)

df_mouse_long <- df_mouse %>%
  select(ID_to_plot, celltype_level_1, percent) %>%
  pivot_longer(cols = percent, names_to = "name", values_to = "value")

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 13), 
        legend.position = "right")

ggplot(df_mouse_long, aes(x = celltype_level_1, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")


group_map <- Ktx_data_mouse@meta.data %>% distinct(ID_to_plot, group)
df_mouse <- df_mouse %>% left_join(group_map, by = "ID_to_plot")

# sufficient to only test Leuko
leuko_ttest <- df_mouse %>%
  filter(celltype_level_1 == "Leuko") %>%
  group_by(celltype_level_1) %>%  
  t_test(percent ~ group, var.equal = FALSE) %>%  
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(celltype_level_1, statistic, df, p, p.adj, p.adj.signif, 
         group1, group2, n1, n2) %>%
  ungroup()

print(leuko_ttest)


