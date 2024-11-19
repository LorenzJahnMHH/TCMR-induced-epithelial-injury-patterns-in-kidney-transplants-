library(Seurat)
library(harmony)

setwd("D:/R_Ordner/Ktx Daten")
Ktx_data_human <- readRDS("Ktx_human_int_harmony.rds") 

# PT subclustering human 
PT_cells <- subset(Ktx_data_human, subset = broad_celltype == c("PT"))

obj.list <- SplitObject(PT_cells, split.by = "sample")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors_human <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

PT_human_KTX_subclustering <- IntegrateData(anchorset = seurat.anchors_human)


PT_human_KTX_subclustering <- ScaleData(PT_human_KTX_subclustering, verbose = FALSE)
PT_human_KTX_subclustering <- RunPCA(PT_human_KTX_subclustering, npcs = 30, verbose = FALSE)
PT_human_KTX_subclustering <- FindNeighbors(PT_human_KTX_subclustering, reduction = "pca", dims = 1:10)
PT_human_KTX_subclustering <- FindClusters(PT_human_KTX_subclustering, resolution = 0.5)
PT_human_KTX_subclustering <- RunUMAP(PT_human_KTX_subclustering, reduction = "pca", dims = 1:10)
DimPlot(PT_human_KTX_subclustering, label = T)

# TAL subclustering human 
TAL_cells_human <- subset(Ktx_human_int_harmony, subset = broad_celltype == c("TAL"))

obj.list <- SplitObject(TAL_cells_human, split.by = "sample")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

TAL_human_KTX_subclustering <- IntegrateData(anchorset = seurat.anchors, k.weight = 40)


TAL_human_KTX_subclustering <- ScaleData(TAL_human_KTX_subclustering, verbose = FALSE)
TAL_human_KTX_subclustering <- RunPCA(TAL_human_KTX_subclustering, npcs = 30, verbose = FALSE)
TAL_human_KTX_subclustering <- FindNeighbors(TAL_human_KTX_subclustering, reduction = "pca", dims = 1:10)
TAL_human_KTX_subclustering <- FindClusters(TAL_human_KTX_subclustering, resolution = 0.5)
TAL_human_KTX_subclustering <- RunUMAP(TAL_human_KTX_subclustering, reduction = "pca", dims = 1:10)
DefaultAssay(TAL_human_KTX_subclustering) <- "RNA"
DimPlot(TAL_human_KTX_subclustering, label = T)

