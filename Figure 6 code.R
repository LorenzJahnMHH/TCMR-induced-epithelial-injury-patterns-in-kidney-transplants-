# this code presents the approach to determine the gene sets used in the analysis of the bulk transcriptomics samples

# Function to compute local densities block-wise
compute_blockwise <- function(expression_data, neighbor_matrix, block_size = 1000) {
  results <- list()
  num_genes <- nrow(expression_data)
  for (start in seq(1, num_genes, by = block_size)) {
    end <- min(start + block_size - 1, num_genes)
    block <- expression_data[start:end, , drop = FALSE]
    results[[length(results) + 1]] <- block %*% neighbor_matrix
  }
  do.call(rbind, results)
}

genGeneSetCount <- function(seurat_object, selected_genes) {
  # Extract raw counts only for the selected genes
  raw_counts <- GetAssayData(seurat_object, slot = "counts")[selected_genes, , drop = FALSE]
  
  # Compute column sums (library sizes) once
  library_sizes <- colSums(GetAssayData(seurat_object, slot = "counts"))
  
  # Calculate CPM only for selected genes
  cpm <- sweep(raw_counts, 2, library_sizes, FUN = "/") * 1e6
  
  # Average CPM for selected genes
  average_cpm <- colMeans(cpm)
  
  # Add to meta.data
  seurat_object@meta.data$average_cpm <- average_cpm
  return(seurat_object)
}

library(FNN) # For nearest neighbor calculations
library(Seurat)
library(Matrix)
library(readxl)
library(grid)

ktx = readRDS("Ktx_data_human.rds") # full human snRNA-seq data
pt=readRDS("PT_human_KTX_subclustering.rds")
tal=readRDS("TAL_human_KTX_subclustering.rds")
gene.df.list = list()
ret = data.frame(matrix(NA, nrow = 100, ncol = 4))
colnames(ret) = c("PT_Injury_h1","PT_Injury_h2","TAL_Injury_h1","TAL_Injury_h2")
ktx$cois <- ifelse(ktx$celltype_level_2 %in% colnames(ret), ktx$celltype_level_2, "Other")

# Step 1: Extract UMAP coordinates and expression matrix
umap_coords <- Embeddings(ktx, "umap")
pairs=list(c("TAL_Injury_h1","TAL_Injury_h2"), c("TAL_Injury_h2", "TAL_Injury_h1"),
        c("PT_Injury_h1","PT_Injury_h2"), c("PT_Injury_h2","PT_Injury_h1"))

so.coff=1.25
socoi.coff=1
gene.list = list()
ret = data.frame(matrix(NA, nrow = 100, ncol = 4))
colnames(ret) = c("PT_Injury_h1","PT_Injury_h2","TAL_Injury_h1","TAL_Injury_h2")
tag = paste0(so.coff,".",socoi.coff)
# Step 2: Identify small and other cells
pdf(paste0("markers.",
           tag,".pdf"), onefile = TRUE)
for(p in pairs){
  print(p)
  coi = p[1]
  coi.oth=p[2]
  small_cells <- rownames(ktx@meta.data[ktx@meta.data$celltype_level_2 == coi, ])
  other_cells <- setdiff(rownames(umap_coords), small_cells)
  
  #compute markers--
  Idents(ktx) <- ifelse(Cells(ktx) %in% small_cells, "Small_Pop", "Other")
  markers <- FindMarkers(
    ktx, 
    ident.1 = "Small_Pop", 
    ident.2 = "Other", 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    only.pos=T
  )
  genes.tmp = rownames(markers[markers$p_val_adj<0.05,])
  length(genes.tmp)
  expression_data <- GetAssayData(ktx, slot = "data", )[genes.tmp,]
  
  av.small.cells <- GetAssayData(ktx, slot = "data")[genes.tmp, small_cells, drop = FALSE]
  av.small.cells <- rowMeans(av.small.cells)
  
  cells.tmp = rownames(ktx@meta.data[ktx@meta.data$celltype_level_2==coi.oth,])
  av.coi.oth <- GetAssayData(ktx, slot = "data")[genes.tmp, cells.tmp, drop = FALSE]
  av.coi.oth <- rowMeans(av.coi.oth)
  
  # Step 3: Exclude neighboring cells
  k = min(100,length(small_cells)-1)
  nn_to_small_pop <- get.knn(umap_coords[small_cells, ], k = k) # Adjust k as needed
  neighboring_cells <- unique(as.vector(nn_to_small_pop$nn.index))
  filtered_other_cells <- setdiff(other_cells, neighboring_cells)
  
  # Step 4: Calculate gene densities
  filtered_coords <- umap_coords[filtered_other_cells, ]
  expression_data_filtered <- as.matrix(expression_data[, filtered_other_cells]) # Convert sparse -> dense
  
  # Efficient neighbor matrix
  nn <- get.knn(filtered_coords, k = k)
  neighbor_indices <- nn$nn.index
  
  # Create neighbor matrix (binary adjacency matrix)
  neighbor_matrix_sparse <- sparseMatrix(
    i = rep(1:nrow(neighbor_indices), each = ncol(neighbor_indices)),
    j = as.vector(neighbor_indices),
    x = 1,
    dims = c(nrow(filtered_coords), nrow(filtered_coords))
  )
  
  # Normalize neighbor_matrix to get mean (divide by k)
  neighbor_matrix_sparse <- neighbor_matrix_sparse / k
  
  # Step 2: Convert expression_data_filtered to a sparse matrix
  expression_data_sparse <- Matrix(expression_data_filtered, sparse = TRUE)
  rm(expression_data_filtered)
  gc()
  
  # Run block-wise computation
  local_densities <- compute_blockwise(expression_data_sparse, neighbor_matrix_sparse)
  
  # Aggregate density per gene
  #gene_densities <- rowMeans(local_densities)
  gene_densities <- apply(local_densities,1, max)
  gene_densities = sort(gene_densities, decreasing = T)
  head(gene_densities)
  tail(gene_densities)
  
  gene.df <- data.frame(
    other = gene_densities[names(gene_densities)],
    small = av.small.cells[names(gene_densities)],
    coi.oth = av.coi.oth,
    ratio.small.other = av.small.cells[names(gene_densities)]/gene_densities[names(gene_densities)],
    ratio.small.coi.oth = av.small.cells[names(gene_densities)]/av.coi.oth[names(gene_densities)]
  )
  
  gene.df = gene.df[order(gene.df$ratio.small.other, decreasing = T),]
  genes.tmp = rownames(gene.df[gene.df$ratio.small.other>so.coff&gene.df$ratio.small.coi.oth>socoi.coff,])
  length(genes.tmp)
  genes.tmp = genes.tmp[1:min(100,length(genes.tmp))]
  gene.list[[coi]] = genes.tmp
  ret[[coi]] = c(genes.tmp, rep(NA, nrow(ret) - length(genes.tmp)))
  
  ktx = genGeneSetCount(ktx, genes.tmp)
  ktx$custom_color <- ifelse(Cells(ktx) %in% small_cells, "Small Cells", "Other")
  
  grid.newpage()
  grid.text(paste("Plots for:", coi, "#genes:", length(genes.tmp)), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  print(FeaturePlot(ktx, features = "average_cpm"))
  print(VlnPlot(ktx, features = "average_cpm", group.by = "cois", sort = TRUE))
  print(DimPlot(ktx, group.by = "custom_color", cols = c("grey", "red")))
}
dev.off()

#if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
library(writexl)
write_xlsx(ret, paste0("markers",tag,".xlsx"))


