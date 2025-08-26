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

# Figure 4D

#spatial proximity calculation-----
library(RANN)
library(ggsignif)

sams = c("bl6.bc.3", "bl6.bc.1", "bc.bl6.1", "bl6.bc.2", "bc.bl6.3")
clusts=c("PT_Injury_m1","PT_Injury_m2","PT_Injury_m3","PT_Injury_m4","TAL_Injury_m1",
           "TAL_Injury_m2", "TAL_Injury_m3","prolifTAL","prolifPT","random")

ctois = c("Cd8+ T cells","Cd4+ T cells","Myeloid","NK","B_cells")

mat = as.data.frame(matrix(-1,ncol=4, nrow=length(sams)*length(clusts)*length(ctois)))
colnames(mat) = c("sam","inj.clust","leuk.clust","perc.dir.nb")
mat$sam = rep(sams, each=length(clusts)*length(ctois))
mat$inj.clust = rep(clusts, length(sams)*length(ctois))
mat$leuk.clust = rep(rep(ctois, each=length(clusts)),length(sams)) 
rownames(mat) = paste(mat$sam,mat$inj.clust,mat$leuk.clust,sep = "_")
mat.list = list()

wfo = "/my_dir/"
set.seed(0)
radius = 25
for(sam in sams){
  print(sam)
  mat.list[[sam]] = list()
  
  tmp = as.data.frame(read.table(paste0(wfo,sam,".cells.csv"), 
                                 sep = ",")) # output from xenium ranger
  colnames(tmp) = tmp[1,]
  tmp = tmp[-1,]
  rownames(tmp) = tmp$cell_id
  cells.tmp = intersect(rownames(tmp), spat@meta.data$orig.cellname)
  tmp = tmp[cells.tmp,]
  
  #transferring celltype labels to tmp
  cells.tmp = rownames(spat@meta.data[spat@meta.data$sample==sam,])
  tmp2 = as.data.frame(spat@meta.data[cells.tmp,c("celltype2","orig.cellname")])
  cells.tmp = spat@meta.data[cells.tmp,]$orig.cellname
  cells.tmp = intersect(cells.tmp, rownames(tmp))
  tmp2 = tmp2[tmp2$orig.cellname %in% cells.tmp,]
  tmp2$cellname = rownames(tmp2)
  rownames(tmp2) = tmp2$orig.cellname
  
  tmp = tmp[cells.tmp,]
  tmp$ct = tmp2[rownames(tmp),"celltype2"]
  tmp$x_centroid = as.numeric(tmp$x_centroid)
  tmp$y_centroid = as.numeric(tmp$y_centroid)
  #choose random cells
  cells.tmp = sample(rownames(tmp[!(tmp$ct %in% ctois), ]),1000, replace=F)
  tmp[cells.tmp,]$ct = "random"
  rm(tmp2)
  #gc()
  
  for(ctoi in ctois){
    #print(ctoi)
    tmp[[ctoi]] = -1
    cells.tmp = rownames(tmp[tmp$ct==ctoi,])  
    tmp3 = tmp[cells.tmp,c("x_centroid","y_centroid")]
    tmp3 = as.matrix(tmp3)
    
    cells.tmp = rownames(tmp[tmp$ct %in% clusts,])
    tmp4 = tmp[cells.tmp,c("x_centroid","y_centroid")]
    tmp4 = as.matrix(tmp4)
    
    nn_result <- nn2(data = tmp3, query = tmp4, k = 1)
    
    # Minimum distances
    min_distances <- nn_result$nn.dists
    
    tmp[cells.tmp,ctoi] = min_distances
    
    for(clust in clusts){
      tmp3 = tmp[tmp$ct==clust,][[ctoi]]
      mat[paste(sam,clust,ctoi,sep = "_"),"perc.dir.nb"] = length(tmp3[tmp3<radius])/length(tmp3)
    }
  }
  mat.list[[sam]] = tmp 
}

#renaming
clusts.transl=c("PT_Injury_m1"="PT Injury m1","PT_Injury_m2"="PT Injury m2",
                "PT_Injury_m3"="PT Injury m3","PT_Injury_m4"="PT Injury m4",
                "TAL_Injury_m1"="TAL Injury m1","TAL_Injury_m2"="TAL Injury m2", 
                "TAL_Injury_m3"="TAL Injury m3","prolifTAL"="TAL Prolif",
                "prolifPT"="PT Prolif","random"="random cells")
ctois = c("Cd8+ T cells","Cd4+ T cells","Myeloid","NK","B_cells")
mat$inj.clust = as.vector(clusts.transl[mat$inj.clust])

clusts=c("PT Prolif","PT Injury m1","PT Injury m2","PT Injury m3","PT Injury m4","TAL Prolif",
         "TAL Injury m1","TAL Injury m2","TAL Injury m3","random cells")
mat$inj.clust = factor(mat3$inj.clust, levels = clusts)

plot.list = list()
for(ctoi in ctois){
  mat3 = mat[mat$leuk.clust==ctoi & mat$inj.clust %in% clusts,]
  mat3$perc.dir.nb = as.numeric(mat3$perc.dir.nb)*100
  
  # Ensure there are no unused levels in the 'sam' variable
  
  # Plot
  plot.list[[ctoi]] = ggplot(mat3, aes(x = inj.clust, y = perc.dir.nb)) + 
    geom_boxplot(fill = NA, outlier.shape = NA) + 
    geom_jitter(aes(shape = sam), color = "black", fill = "black", size = 1.2, alpha = 0.9) + 
    scale_shape_manual(values = c("bc.bl6.1" = 24, "bc.bl6.3" = 24, "bl6.bc.1" = 25, "bl6.bc.2" = 25, "bl6.bc.3" = 25, "bc.bc" = 16, "bl6.bl6" = 15)) + 
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
      axis.title.y = element_text(size = 9),
      panel.background = element_blank(),  # Remove background color
      plot.background = element_blank()
    ) + 
    ggtitle("") + 
    xlab("") + 
    ylab("% cells with direct neighbour") +
    expand_limits(y = 0) +
    ggtitle(ctoi)
  
}

library(gridExtra)

#plot 1
ctois = c("Cd8+ T cells","Cd4+ T cells","Myeloid")
a4_width <- 8.27
a4_height <- 11.69 / 5  
plot.list <- lapply(seq_along(plot.list), function(i) {
  if (i == 1) {
    # Keep the y-axis label and decrease its font size for the leftmost plot
    plot.list[[i]] + theme(axis.title.y = element_text(size = 8))
  } else {
    # Remove the y-axis label for other plots
    plot.list[[i]] + theme(axis.title.y = element_blank())
  }
})
pdf("/Users/chhinze/Dropbox/Projects/dr.arbeit.jahn/paper.spatial/figures/spat.dist.leuk.1.pdf", width = a4_width, height = a4_height)
grid.arrange(grobs = plot.list, ncol = length(plot.list))
dev.off()




                               

