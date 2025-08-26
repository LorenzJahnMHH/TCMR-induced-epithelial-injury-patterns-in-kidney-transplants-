# R v4.4.3
library(reticulate) # v1.42.0
library(scCustomize) # v3.0.1
library(dplyr) # v1.1.4
library(purrr) # v1.0.4
library(tibble) # v3.2.1

# Seurat v5.3.0
# SeuratObject v5.1.0

# Python v3.13.3
# AnnData v0.11.4

use_condaenv("seurat5")

data <- readRDS("./spat_reseg5.rds")

data[["BlankCodeword"]] <- NULL
data[["ControlCodeword"]] <- NULL
data[["ControlProbe"]] <- NULL

transl <- c("bc.bc" = "control", "bc.bl6.1" = "rejection", "bc.bl6.3" = "rejection", "bl6.bc.1" = "rejection", "bl6.bc.2" = "rejection", "bl6.bc.3" = "rejection", "bl6.bl6" = "control")
data$broad.group <- as.vector(transl[data$sample])

transl <- c("bc.bc" = "control.mild.rejection", "bc.bl6.1" = "mild.rejection", "bc.bl6.3" = "mild.rejection", "bl6.bc.1" = "strong.rejection", "bl6.bc.2" = "strong.rejection", "bl6.bc.3" = "strong.rejection", "bl6.bl6" = "control.strong.rejection")
data$finer.group <- as.vector(transl[data$sample])

coords <- map(Images(data), ~ GetTissueCoordinates(data[[.]])) |>
  bind_rows() |>
  column_to_rownames("cell")

data[[c("x", "y")]] <- coords[colnames(data), c("x", "y")]

as.anndata(x = data, file_path = ".", file_name = "spat_reseg5.h5ad")
