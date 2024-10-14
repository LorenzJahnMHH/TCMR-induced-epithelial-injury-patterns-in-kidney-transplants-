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

# Leuko subclustering human 
Leukos_obj <- subset(Ktx_human_int_harmony, subset = broad_celltype == c("Leuko")) # bzw. Leukocyte_humanKTX_subclustering
obj.list <- SplitObject(Leukos_obj, split.by = "group")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Integrated_Leukos_human <- IntegrateData(anchorset = seurat.anchors)



Integrated_Leukos_human <- ScaleData(Integrated_Leukos_human, verbose = FALSE)
Integrated_Leukos_human <- RunPCA(Integrated_Leukos_human, npcs = 30, verbose = FALSE)
Integrated_Leukos_human <- FindNeighbors(Integrated_Leukos_human, reduction = "pca", dims = 1:15)
Integrated_Leukos_human <- FindClusters(Integrated_Leukos_human, resolution = 0.5)
Integrated_Leukos_human <- RunUMAP(Integrated_Leukos_human, reduction = "pca", dims = 1:15)
DefaultAssay(Integrated_Leukos_human) <- "RNA"
DimPlot(Integrated_Leukos_human, group.by = 'seurat_clusters', label = T) 

Integrated_Leukos_human_harmony2 <- harmony::RunHarmony(Integrated_Leukos_human, group.by.vars = c("sample"), lambda = 2, tau = 1000, theta = 1, assay.use = "integrated",
                                                        kmeans_init_nstart = 60, kmeans_init_iter_max=2000,
                                                        max.iter.harmony = 30)


Integrated_Leukos_human_harmony2 <- ScaleData(Integrated_Leukos_human_harmony2, verbose = FALSE)
Integrated_Leukos_human_harmony2 <- FindNeighbors(Integrated_Leukos_human_harmony2, reduction = "harmony", dims = 1:15)
Integrated_Leukos_human_harmony2 <- FindClusters(Integrated_Leukos_human_harmony2, resolution = 0.5)


library(harmony)
DefaultAssay(Leukos_obj) <- "RNA"  # Make sure RNA is the active assay if Harmony needs it as input

# Pre-process Leukos_obj as needed
Leukos_obj <- NormalizeData(Leukos_obj)
Leukos_obj <- FindVariableFeatures(Leukos_obj)
Leukos_obj <- ScaleData(Leukos_obj, verbose = FALSE)
Leukos_obj <- RunPCA(Leukos_obj, npcs = 30, verbose = FALSE)

Leukos_obj <- RunHarmony(
  Leukos_obj,
  group.by.vars = "sample",  # or "group" if that's more accurate for your data
  reduction = "pca",
  dims = 1:15,  # based on PCA dimensions to include
  lambda = 2, 
  tau = 1000, 
  theta = 1
)

Leukos_obj <- ScaleData(Leukos_obj, verbose = FALSE)
Leukos_obj <- FindNeighbors(Leukos_obj, reduction = "harmony", dims = 1:15)
Leukos_obj <- FindClusters(Leukos_obj, resolution = 0.5)
Leukos_obj <- RunUMAP(Leukos_obj, reduction = "harmony", dims = 1:15)

# Visualize the final integrated data
DimPlot(Leukos_obj, group.by = "seurat_clusters", label = TRUE)
