# python = 3.10.14
# anndata = 0.10.7
# pandas = 2.2.2

import anndata as ad
import pandas as pd

filepath = "spat_reseg5.h5ad"
metadata = "meta.data.txt"
pt_meta = "pt.meta.txt"
tal_meta = "tal.meta.txt"

celltype_annotation = pd.read_table(metadata, index_col=0)

# level_2 and level_3 are identical for PT and TAL
new_annotation = pd.concat(
    [pd.read_table(f, index_col=0) for f in [pt_meta, tal_meta]]
).assign(celltype_level_2=lambda df: df["celltype_level_3"])

# update with corrected annotation
celltype_annotation.update(new_annotation[["celltype_level_2", "celltype_level_3"]])

# prepare and only keep necessary information
adata = ad.read_h5ad(filepath)

del adata.layers

adata.var = adata.var[[]]

adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy().astype("float32")

adata.obs = adata.obs[
    ["sample", "broad.group", "finer.group", "celltype_level_1"]
].join(celltype_annotation)

adata.write_h5ad(filepath)
