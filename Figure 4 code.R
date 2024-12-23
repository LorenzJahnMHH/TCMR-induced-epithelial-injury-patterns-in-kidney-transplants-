# Figure 4
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

setwd("D:/R_Ordner/Ktx Daten")
Ktx_Leuko_mouse_subclustering <- readRDS("D:/R_Ordner/KTx Daten/Ktx_Leukos_mouse_subclustering_harmony.rds")

# Figure 4A UMAP

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

# Figure 4A Heatmap
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




                               
