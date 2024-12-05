setwd("D:/R_Ordner/KTx Daten")
Ktx_data <- readRDS("Ktx_data.rds")

library(Seurat)
library(dplyr)

#DCT
Prolif_subset <- subset(
  Ktx_data,
  subset = celltype_level_1 == "Prolif" & UMAP_1 > 6.5 & UMAP_1 < 10 & UMAP_2 < -10.25) # keep DCT cells from Prolif cluster
DCT_cells <- subset(Ktx_data, subset = celltype_level_1 == "DCT")

Combined_subset <- merge(DCT_cells, y = Prolif_subset)


obj.list <- SplitObject(Combined_subset, split.by = "ID")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_DCT_mouse <- IntegrateData(anchorset = seurat.anchors)
Integrated_DCT_mouse <- ScaleData(Integrated_DCT_mouse, verbose = FALSE)
Integrated_DCT_mouse <- RunPCA(Integrated_DCT_mouse, npcs = 30, verbose = FALSE)
Integrated_DCT_mouse <- FindNeighbors(Integrated_DCT_mouse, reduction = "pca", dims = 1:10)
Integrated_DCT_mouse <- FindClusters(Integrated_DCT_mouse, resolution = 0.5)
Integrated_DCT_mouse <- RunUMAP(Integrated_DCT_mouse, reduction = "pca", dims = 1:10)

#CNT 
CNT_cells <- subset(Ktx_data, subset = celltype_level_1 == c("CNT"))
obj.list <- SplitObject(CNT_cells, split.by = "ID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_CNT_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 79)


Integrated_CNT_mouse <- ScaleData(Integrated_CNT_mouse, verbose = FALSE)
Integrated_CNT_mouse <- RunPCA(Integrated_CNT_mouse, npcs = 30, verbose = FALSE)
Integrated_CNT_mouse <- FindNeighbors(Integrated_CNT_mouse, reduction = "pca", dims = 1:10)
Integrated_CNT_mouse <- FindClusters(Integrated_CNT_mouse, resolution = 0.5)
Integrated_CNT_mouse <- RunUMAP(Integrated_CNT_mouse, reduction = "pca", dims = 1:10)
# exclude doublet Cluster 5 
Integrated_CNT_mouse <- subset(Integrated_CNT_mouse, subset = seurat_clusters == "5", invert = T)

obj.list <- SplitObject(Integrated_CNT_mouse, split.by = "ID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_CNT_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 77)


Integrated_CNT_mouse <- ScaleData(Integrated_CNT_mouse, verbose = FALSE)
Integrated_CNT_mouse <- RunPCA(Integrated_CNT_mouse, npcs = 30, verbose = FALSE)
Integrated_CNT_mouse <- FindNeighbors(Integrated_CNT_mouse, reduction = "pca", dims = 1:30)
Integrated_CNT_mouse <- FindClusters(Integrated_CNT_mouse, resolution = 0.5)
Integrated_CNT_mouse <- RunUMAP(Integrated_CNT_mouse, reduction = "pca", dims = 1:30)

## EC
EC_cells <- subset(Ktx_data, subset = celltype_level_1 == c("EC"))
obj.list <- SplitObject(EC_cells, split.by = "ID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_EC_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 33)


Integrated_EC_mouse <- ScaleData(Integrated_EC_mouse, verbose = FALSE)
Integrated_EC_mouse <- RunPCA(Integrated_EC_mouse, npcs = 30, verbose = FALSE)
Integrated_EC_mouse <- FindNeighbors(Integrated_EC_mouse, reduction = "pca", dims = 1:30)
Integrated_EC_mouse <- FindClusters(Integrated_EC_mouse, resolution = 0.5)
Integrated_EC_mouse <- RunUMAP(Integrated_EC_mouse, reduction = "pca", dims = 1:30)
# exclude doublet Cluster 3 + 7 (isolated doublets) 

#tL
tL_cells <- subset(Ktx_data, subset = celltype_level_1 == c("tL"))
obj.list <- SplitObject(tL_cells, split.by = "group")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_tL_mouse <- IntegrateData(anchorset = seurat.anchors)


Integrated_tL_mouse <- ScaleData(Integrated_tL_mouse, verbose = FALSE)
Integrated_tL_mouse <- RunPCA(Integrated_tL_mouse, npcs = 30, verbose = FALSE)
Integrated_tL_mouse <- FindNeighbors(Integrated_tL_mouse, reduction = "pca", dims = 1:30)
Integrated_tL_mouse <- FindClusters(Integrated_tL_mouse, resolution = 0.5)
Integrated_tL_mouse <- RunUMAP(Integrated_tL_mouse, reduction = "pca", dims = 1:30)

# exclude doublet Cluster 5 + 6 
Integrated_tL_mouse <- subset(Integrated_tL_mouse, subset = seurat_clusters %in%  c("5", "6"), invert = T)

obj.list <- SplitObject(Integrated_tL_mouse, split.by = "group")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_tL_mouse <- IntegrateData(anchorset = seurat.anchors)


Integrated_tL_mouse <- ScaleData(Integrated_tL_mouse, verbose = FALSE)
Integrated_tL_mouse <- RunPCA(Integrated_tL_mouse, npcs = 30, verbose = FALSE)
Integrated_tL_mouse <- FindNeighbors(Integrated_tL_mouse, reduction = "pca", dims = 1:30)
Integrated_tL_mouse <- FindClusters(Integrated_tL_mouse, resolution = 0.5)
Integrated_tL_mouse <- RunUMAP(Integrated_tL_mouse, reduction = "pca", dims = 1:30)

#IntC
IntC_cells <- subset(Ktx_data, subset = celltype_level_1 == c("IntC"))
obj.list <- SplitObject(IntC_cells, split.by = "group") 


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_IntC_mouse <- IntegrateData(anchorset = seurat.anchors)


Integrated_IntC_mouse <- ScaleData(Integrated_IntC_mouse, verbose = FALSE)
Integrated_IntC_mouse <- RunPCA(Integrated_IntC_mouse, npcs = 30, verbose = FALSE)
Integrated_IntC_mouse <- FindNeighbors(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)
Integrated_IntC_mouse <- FindClusters(Integrated_IntC_mouse, resolution = 0.5)
Integrated_IntC_mouse <- RunUMAP(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)

# exclude doublet Cluster 6 + 7 
Integrated_IntC_mouse <- subset(Integrated_IntC_mouse, subset = seurat_clusters %in%  c("6", "7"), invert = T)

obj.list <- SplitObject(Integrated_IntC_mouse, split.by = "group") 


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_IntC_mouse <- IntegrateData(anchorset = seurat.anchors)


Integrated_IntC_mouse <- ScaleData(Integrated_IntC_mouse, verbose = FALSE)
Integrated_IntC_mouse <- RunPCA(Integrated_IntC_mouse, npcs = 30, verbose = FALSE)
Integrated_IntC_mouse <- FindNeighbors(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)
Integrated_IntC_mouse <- FindClusters(Integrated_IntC_mouse, resolution = 0.5)
Integrated_IntC_mouse <- RunUMAP(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)

# exclude doublet Cluster 1 
Integrated_IntC_mouse <- subset(Integrated_IntC_mouse, subset = seurat_clusters == "1", invert = T)

obj.list <- SplitObject(Integrated_IntC_mouse, split.by = "group") 


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_IntC_mouse <- IntegrateData(anchorset = seurat.anchors)


Integrated_IntC_mouse <- ScaleData(Integrated_IntC_mouse, verbose = FALSE)
Integrated_IntC_mouse <- RunPCA(Integrated_IntC_mouse, npcs = 30, verbose = FALSE)
Integrated_IntC_mouse <- FindNeighbors(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)
Integrated_IntC_mouse <- FindClusters(Integrated_IntC_mouse, resolution = 0.5)
Integrated_IntC_mouse <- RunUMAP(Integrated_IntC_mouse, reduction = "pca", dims = 1:30)


#CD-PC
PC_cells <- subset(Ktx_data, subset = celltype_level_1 == c("CD-PC"))
obj.list <- SplitObject(PC_cells, split.by = "groupID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_PC_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 52) 


Integrated_PC_mouse <- ScaleData(Integrated_PC_mouse, verbose = FALSE)
Integrated_PC_mouse <- RunPCA(Integrated_PC_mouse, npcs = 30, verbose = FALSE)
Integrated_PC_mouse <- FindNeighbors(Integrated_PC_mouse, reduction = "pca", dims = 1:10)
Integrated_PC_mouse <- FindClusters(Integrated_PC_mouse, resolution = 0.8)  #resolution was increased to precisely remove doublet clusters
Integrated_PC_mouse <- RunUMAP(Integrated_PC_mouse, reduction = "pca", dims = 1:10)

# exclude doublet Cluster 5 + 6 
Integrated_PC_mouse <- subset(Integrated_PC_mouse, subset = seurat_clusters %in% c("5", "6"), invert = T)

obj.list <- SplitObject(Integrated_PC_mouse, split.by = "groupID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_PC_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 43) 


Integrated_PC_mouse <- ScaleData(Integrated_PC_mouse, verbose = FALSE)
Integrated_PC_mouse <- RunPCA(Integrated_PC_mouse, npcs = 30, verbose = FALSE)
Integrated_PC_mouse <- FindNeighbors(Integrated_PC_mouse, reduction = "pca", dims = 1:10)
Integrated_PC_mouse <- FindClusters(Integrated_PC_mouse, resolution = 0.5) #resolution back to 0.5
Integrated_PC_mouse <- RunUMAP(Integrated_PC_mouse, reduction = "pca", dims = 1:10)

# CD-IC-A CD-IC-B
IC_cells <- subset(Ktx_data, subset = celltype_level_1 %in% c("CD-IC-A", "CD-IC-B"))
obj.list <- SplitObject(IC_cells, split.by = "ID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_IC_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 64) 


Integrated_IC_mouse <- ScaleData(Integrated_IC_mouse, verbose = FALSE)
Integrated_IC_mouse <- RunPCA(Integrated_IC_mouse, npcs = 30, verbose = FALSE)
Integrated_IC_mouse <- FindNeighbors(Integrated_IC_mouse, reduction = "pca", dims = 1:30)
Integrated_IC_mouse <- FindClusters(Integrated_IC_mouse, resolution = 0.5)
Integrated_IC_mouse <- RunUMAP(Integrated_IC_mouse, reduction = "pca", dims = 1:30)

# exclude doublet Cluster 4 + 5 
Integrated_IC_mouse <- subset(Integrated_IC_mouse, subset = seurat_clusters %in% c("4", "5"), invert = T)

obj.list <- SplitObject(Integrated_IC_mouse, split.by = "ID")


for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_IC_mouse <- IntegrateData(anchorset = seurat.anchors, k.weight = 60) 


Integrated_IC_mouse <- ScaleData(Integrated_IC_mouse, verbose = FALSE)
Integrated_IC_mouse <- RunPCA(Integrated_IC_mouse, npcs = 30, verbose = FALSE)
Integrated_IC_mouse <- FindNeighbors(Integrated_IC_mouse, reduction = "pca", dims = 1:30)
Integrated_IC_mouse <- FindClusters(Integrated_IC_mouse, resolution = 0.5)
Integrated_IC_mouse <- RunUMAP(Integrated_IC_mouse, reduction = "pca", dims = 1:30)
DefaultAssay(Integrated_IC_mouse) <- "RNA"
DimPlot(Integrated_IC_mouse, label = T)