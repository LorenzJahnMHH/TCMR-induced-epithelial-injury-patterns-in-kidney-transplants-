# file describing initial clustering of spatial data after xeniumranger and cell type annotation

# initial clustering
library(Seurat)
library(harmony)
options(future.globals.maxSize=1e20)
sl = read.table("sample.list.txt")
ret=list()
for(el in rownames(sl)){
  sam = sl[el,1]
  print(sam)
  fo = sl[el,2]
  ret[[sam]] = LoadXenium(fo)
}

spat = merge(ret[[1]],c(ret[[2]],ret[[3]],ret[[4]],ret[[5]],ret[[6]],ret[[7]]), add.cell.ids=sl[,1], project="spat")
saveRDS(spat, file="spat.rds")
spat = NormalizeData(spat)
spat = FindVariableFeatures(spat)
spat = ScaleData(spat)
spat = RunPCA(spat)
spat@meta.data$orig.ident = unlist(lapply(rownames(spat@meta.data),function(x) strsplit(x, "_")[[1]][1]))
spat = RunHarmony(spat, group.by.vars="orig.ident")

spat <- RunUMAP(spat, reduction = "harmony", dims = 1:30)
spat <- FindNeighbors(spat, reduction = "harmony", dims = 1:30)
spat <- FindClusters(spat, resolution = 0.3)

saveRDS(spat, file="spat.rds")

# major cell type annotation
                                    
spat$celltype = ""
cts = c("0"="non-tub","1"="PT","2"="TL","3"="non-tub","4"="TAL",
        "5"="PT","6"="non-tub","7"="non-tub","8"="CD-PC","9"="DCT",
        "10"="non-tub","11"="Podo","12"="non-tub","13"="Uro","14"="non-tub","15"="non-tub")
# this annotation was determined by investigating marker gene expression and distribution of cells on the kidney slide
spat$celltype = as.vector(cts[as.vector(spat$Xenium_snn_res.0.3)])

spat@meta.data$celltype2 = spat@meta.data$celltype
# non-tub sub annotation (e.g. EC, interstitial cells, leukocytes) using RCTD
cells.tmp = rownames(spat@meta.data[spat@meta.data$Xenium_snn_res.0.3 %in% c("0","3","6","7","10","12","14"),])

library(spacexr)
ktx <- readRDS("Ktx_data.rds") # this is the mouse snRNA-seq data
ktx <- UpdateSeuratObject(ktx)
Idents(ktx) <- "broad_celltype"

ktx = subset(ktx, idents = c("Endothelial Cells", "Interstitial Cells", "Immune Cells"))
table(Idents(ktx))
ktx$broad_celltype = factor(ktx$broad_celltype)

counts <- GetAssayData(ktx, assay = "RNA", slot = "counts")
cluster <- as.factor(ktx$broad_celltype)
names(cluster) <- colnames(ktx)
nUMI <- ktx$nCount_RNA
names(nUMI) <- colnames(ktx)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

sams = unique(spat$orig.ident)
Idents(spat) = "orig.ident"

ret = data.frame(spot_class=character(), first_type=character(), second_type=character(), first_class=character(), second_class=character(), spot_class=character(), first_type=character(), second_type=character(), first_class=logical(), second_class=logical(), min_score=double(), singlet_score=double(), conv_all=logical(), conv_doublet=logical(), stringsAsFactors=F) 

for(sam in sams){
print(sam)
tmp2 = subset(spat, idents = sam)
query.counts <- GetAssayData(tmp2, assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(tmp2, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 6, UMI_min=0)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
ret = rbind(ret, annotations.df)

results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
}

write.table(ret,file="/my_dir/non-tub.annot.tsv" ,col.names=T, row.names=T, quote=F, sep="\t")

# non-tub.annot.tsv and similar following files are the output of RCTD (robust cell type decomposition).
tmp = read.table(file = "non-tub.annot.tsv", 
                 header = T,
                 sep = "\t")
pie(table(tmp$spot_class))

cells.tmp = rownames(tmp[tmp$spot_class=="singlet",])
spat@meta.data[cells.tmp,]$celltype2 = tmp[cells.tmp,]$first_type
table(spat$celltype2)
cells.tmp = rownames(spat@meta.data[spat@meta.data$celltype2=="non-tub",])
spat@meta.data[cells.tmp,]$celltype2 = "non-tub_na"

# cd-ic, prolif, pec (ct.annot.tsv uses the same approach for non-tub annotation using all level 1 cell types from mouse snRNA-seq and spat.rds to identify the metioned cell types)
tmp = read.table("ct.annot.tsv", header = T, 
                 sep = "\t")
pie(table(tmp$spot_class))

cells.tmp = rownames(tmp[tmp$spot_class=="singlet"&tmp$first_type %in% c("Prolif","CD-IC-A",
                                                                         "CD-IC-B","PEC"),])
pie(table(tmp[cells.tmp,]$first_type))

table(spat@meta.data[cells.tmp,]$celltype2)
spat@meta.data[cells.tmp,]$celltype2 = tmp[cells.tmp,]$first_type
table(spat$celltype2)

#glom - all celltypes annotated as podocyte in spat.rds are re-annotated using mouse snRNA-seq data from podocytes, EC, PEC and IntC 
tmp = read.table("/Users/chhinze/Dropbox/Projects/dr.arbeit.jahn/analysis/subct.annot/glom.annot.tsv", 
                 header = T, 
                 sep = "\t")
pie(table(tmp$spot_class))
pie(table(tmp$first_type))

cells.tmp = rownames(tmp[tmp$spot_class=="singlet",])
spat@meta.data[cells.tmp,]$celltype2 = tmp[cells.tmp,]$first_type
table(spat$celltype2)

cells.tmp = rownames(spat@meta.data[spat@meta.data$celltype2=="Interstitial Cells",])
spat@meta.data[cells.tmp,]$celltype2 = "IntC"
table(spat$celltype2)

cells.tmp = rownames(spat@meta.data[spat@meta.data$celltype2=="Endothelial Cells",])
spat@meta.data[cells.tmp,]$celltype2 = "EC"
table(spat$celltype2)

# further sub annotation of PT, TAL and leukocyte sub cell types in analogous manner using PT, TAL cells or leukocytes from spat and the respective cells with sub cell stae/type annotation from mouse snRNA-seq                                          
