library(Seurat)
library(harmony)

# Clusters with low nCount_RNA and nFeature_RNA were removed from raw data
# Removal of doublet clusters in subclusterings with at least 2 other major celltype canonical marker
# Find the included cells with celltype levels in the metadata sheet available on GEO

obj.list <- SplitObject(Human.samples, split.by = "group")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Ktx_data_human <- IntegrateData(anchorset = seurat.anchors)


Ktx_data_human <- ScaleData(Ktx_data_human, verbose = FALSE)
Ktx_data_human <- RunPCA(Ktx_data_human, npcs = 30, verbose = FALSE)

Ktx_data_human <- harmony::RunHarmony(Ktx_data_human, group.by.vars = c("group"), lambda = 2, tau = 1000, theta = 1, assay.use = "integrated",
                                                    kmeans_init_nstart = 60, kmeans_init_iter_max=2000,
                                                    max.iter.harmony = 30)

Ktx_data_human <- ScaleData(Ktx_data_human, verbose = FALSE)
Ktx_data_human <- FindNeighbors(Ktx_data_human, reduction = "harmony", dims = 1:15)
Ktx_data_human <- FindClusters(Ktx_data_human, resolution = 0.5)
Ktx_data_human <- RunUMAP(Ktx_data_human, reduction = "harmony", dims = 1:15)

# PT subclustering human 
PT_cells <- subset(Ktx_data_human, subset = celltype_level_1 == c("PT"))

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
TAL_cells_human <- subset(Ktx_data_human, subset = celltype_level_1 == c("TAL"))

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
