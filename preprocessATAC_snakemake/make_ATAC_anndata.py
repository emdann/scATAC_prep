import pandas as pd 
import scanpy as sc

bc = pd.read_table(input["bcs"], header=None)
feat = pd.read_table(input["peaks"], header=None)
adata = sc.read_mtx(input["cnts"]).T
adata.obs_names = bc[0]
adata.var_names = feat[0]

## Add peak annotation to var
peak_anno_df = pd.read_csv(input["peak_anno"], index_col=0)
peak_anno_df.index = peak_anno_df["peak_id"]
peak_anno_df.drop("peak_id",1, inplace=True)
adata.var = pd.concat([adata.var, peak_anno_df], 1)

## Save binary data to layers
adata.layers["binary_raw"] = adata.X
adata.layers["binary_raw"][adata.layers["binary_raw"] > 1] = 1

## Filter peaks 
var_qc = sc.pp.calculate_qc_metrics(adata, layer="binary_raw")[1]
adata.var = pd.concat([adata.var, var_qc], 1)
adata = adata[:,adata.var.total_counts > params['k']] # Accessible in at least k cells
adata = adata[:, adata.var.ENCODE_blacklist==0] # Not in ENCODE blacklist
adata = adata[:, adata.var.peak_width > params['width']] # Filter by width (remove peaks at the lowest end)

## Save count matrix for cisTopic
ct_df = pd.DataFrame(adata.layers["binary_raw"].toarray().T)
ct_df.columns = adata.obs_names
ct_df.index = adata.var_names
ct_df.to_csv(output["cnt_tsv"], sep="\t")

## Write output anndata
adata.write_h5ad(output["h5ad"])