library(Seurat)
library(harmony)

# data read in and clustering for cell type annotation
# read in data and create Seurat object for each sample using:
# CreateSeuratObject with nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA < 20000 & percent.mt < 5
# merge all samples in Ktx_data 

Ktx_data <- SplitObject(Ktx_data, split.by = "orig.ident")

Ktx_data <- lapply(X = Ktx_data, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
Ktx_data <- FindIntegrationAnchors(object.list = Ktx_data, dims = 1:30, k.anchor = 5, verbose = TRUE)
Ktx_data <- IntegrateData(anchorset = Ktx_data, dims = 1:30, verbose = TRUE)

DefaultAssay(Ktx_data) <- "integrated"

Ktx_data <- ScaleData(Ktx_data, verbose = TRUE)
Ktx_data <- RunPCA(Ktx_data, npcs = 30, verbose = TRUE)

Ktx_data <- RunUMAP(Ktx_data, reduction = "pca", dims = 1:30, verbose = TRUE)
Ktx_data <- FindNeighbors(Ktx_data, reduction = "pca", dims = 1:30, verbose = TRUE)
Ktx_data <- FindClusters(Ktx_data, resolution = 0.5, verbose = TRUE)
# continue with cell type annotation based on marker gene expression

# PT subclustering mouse 
PT_cells <- subset(Ktx_data, subset = celltype_level_1 == "PT")
obj.list <- SplitObject(PT_cells, split.by = "ID")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Ktx_PT_mouse_subclustering <- IntegrateData(anchorset = seurat.anchors)


Ktx_PT_mouse_subclustering <- ScaleData(Ktx_PT_mouse_subclustering, verbose = FALSE)
Ktx_PT_mouse_subclustering <- RunPCA(Ktx_PT_mouse_subclustering, npcs = 30, verbose = FALSE)
Ktx_PT_mouse_subclustering <- FindNeighbors(Ktx_PT_mouse_subclustering, reduction = "pca", dims = 1:10)
Ktx_PT_mouse_subclustering <- FindClusters(Ktx_PT_mouse_subclustering, resolution = 0.5)
Ktx_PT_mouse_subclustering <- RunUMAP(Ktx_PT_mouse_subclustering, reduction = "pca", dims = 1:10)


# TAL subclustering mouse 
TAL_cells <- subset(Ktx_data, 
  subset = (celltype_level_1 == "TAL") | (celltype_level_2 == "TAL_Prolif")
)


obj.list <- SplitObject(TAL_cells, split.by = "ID")

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Ktx_TAL_mouse_subclustering <- IntegrateData(anchorset = seurat.anchors)


Ktx_TAL_mouse_subclustering <- ScaleData(Ktx_TAL_mouse_subclustering, verbose = FALSE)
Ktx_TAL_mouse_subclustering <- RunPCA(Ktx_TAL_mouse_subclustering, npcs = 30, verbose = FALSE)
Ktx_TAL_mouse_subclustering <- FindNeighbors(Ktx_TAL_mouse_subclustering, reduction = "pca", dims = 1:10)
Ktx_TAL_mouse_subclustering <- FindClusters(Ktx_TAL_mouse_subclustering, resolution = 0.5)
Ktx_TAL_mouse_subclustering <- RunUMAP(Ktx_TAL_mouse_subclustering, reduction = "pca", dims = 1:10)
DefaultAssay(Ktx_TAL_mouse_subclustering) <- "RNA"
DimPlot(Ktx_TAL_mouse_subclustering, label = T)


# Leuko subclustering mouse 
Leuko_cells <- subset(Ktx_data, subset = celltype_level_1 == c("Leuko"))

obj.list <- SplitObject(Leuko_cells, split.by = "donor") # not enough syngeneic cells to split it by ID, correction with harmony after

for(i in 1:length(obj.list)){
  
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

features <- SelectIntegrationFeatures(object.list = obj.list)
seurat.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

Ktx_Leukos_mouse_subclustering <- IntegrateData(anchorset = seurat.anchors)

Ktx_Leuko_mouse_subclustering <- ScaleData(Ktx_Leuko_mouse_subclustering, verbose = FALSE)
Ktx_Leuko_mouse_subclustering <- RunPCA(Ktx_Leuko_mouse_subclustering, npcs = 30, verbose = FALSE)

Ktx_Leuko_mouse_subclustering <- RunHarmony(
  object = Ktx_Leuko_mouse_subclustering, 
  group.by.vars = c("ID"),  
  dims.use = 1:15, 
  assay.use = "integrated",  
  lambda = 2, 
  theta = 1, 
  max.iter.harmony = 30, 
  kmeans_init_nstart = 60, 
  kmeans_init_iter_max = 2000
)

Ktx_Leuko_mouse_subclustering <- RunUMAP(
  object = Ktx_Leuko_mouse_subclustering, 
  reduction = "harmony", 
  dims = 1:15
)

Ktx_Leuko_mouse_subclustering <- FindNeighbors(
  object = Ktx_Leuko_mouse_subclustering, 
  reduction = "harmony", 
  dims = 1:15
)

Ktx_Leuko_mouse_subclustering <- FindClusters(
  object = Ktx_Leuko_mouse_subclustering, 
  graph.name = "integrated_snn",  
  resolution = 0.5
)
# We excluded a doublet with Epithelium markers in the progress 
