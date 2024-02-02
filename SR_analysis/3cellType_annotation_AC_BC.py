import matplotlib.pyplot as plt
import scvelo as scv
import numpy as np
import pandas as pd
import anndata
from os.path import exists
import plotly.express as px
import scanpy as sc
import scvi
import sys
sys.setrecursionlimit(10000)

###integrate short read data for the 4 samples
save_file1 = 'sample1_major.h5ad'
adata1 = sc.read_h5ad(save_file1)
save_file2 = 'sample2_major.h5ad'
adata2 = sc.read_h5ad(save_file2)
save_file3 = 'sample3_major.h5ad'
adata3 = sc.read_h5ad(save_file3)
save_file4 = 'sample4_major.h5ad'
adata4 = sc.read_h5ad(save_file4)

adata1.obs["batch"] = "sample1"
adata2.obs["batch"] = "sample2"
adata3.obs["batch"] = "sample3"
adata4.obs["batch"] = "sample4"

adata = adata1.concatenate(adata2,adata3,adata4,batch_categories=["sample1", "sample2", "sample3", "sample4"])
adata.write_h5ad("SR.h5ad")

###subset bipolar cells
BC_table="cell_anno" #list of all cell barcodes with major cell class annotation 
#barcode	anno
BC = pd.read_csv(BC_table,sep="\t",header=0)
BC.set_index("barcode")
BC.index = BC.barcode
cellname = adata.obs.index.intersection(BC.index)
adata.obs.loc[cellname,"anno"]=BC.loc[cellname,"anno"]
adata_subset = adata[adata.obs['anno'].isin(['BC'])]
#adata_subset = adata[adata.obs['anno'].isin(['AC'])]
adata_subset.write_h5ad("SR_BC_integrated_batch_corrected.h5ad")


###annotate bipolar cell types using scANVI and our in-house meta reference
###for more information for the meta reference, please check https://www.biorxiv.org/content/biorxiv/early/2024/01/28/2024.01.24.577060.full.pdf 
ref="BC_2000.h5ad"  #our in-house metadata reference, or AC_2000.h5ad for amacrine cell type annotation
celltype="BC"  #or AC
ref_label_index="subclass"
queryname="SR"
infile="SR_BC_integrated_batch_corrected.h5ad"
#infile="SR_AC_integrated_batch_corrected.h5ad"

scvi.settings.seed = 0

adata_query=scv.read(infile)
adata_query=adata_query.raw.to_adata()
adata_ref=scv.read(ref)
adata_ref.obs["author_cell_type"]=adata_ref.obs[ref_label_index].astype(str).astype("category")

batch_key = "sample"
labels_key = ""

adata_ref.obs["batch"] = adata_ref.obs.fileid

adata_query.obs["sample_type"] = "query"
adata_query.obs["batch"] = adata_query.obs["sample"]
adata_ref.obs["sample_type"] = "reference"
print("before_concat")


adata_query.var_names=adata_query.var.index.astype("str").astype("category")
adata_ref.var_names=adata_ref.var.index.astype("str").astype("category")
common_genes = adata_ref.var_names.intersection(adata_query.var_names)
adata_ref=adata_ref[:, common_genes]
adata_query=adata_query[:, common_genes]
adata = anndata.concat([adata_ref, adata_query])
print("after_concat")


adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata

batch_key = "batch"

adata.obs[batch_key] = adata.obs[batch_key].astype("str").astype("category")
print(adata.obs[batch_key])
sc.pp.highly_variable_genes(
	adata, n_top_genes=2000, batch_key = batch_key, subset=True
)

adata.obs["celltype_scanvi"]="Unknown"
ss2_mask = adata.obs['sample_type'] == "reference"
adata.obs["celltype_scanvi"][ss2_mask]=adata.obs.author_cell_type[ss2_mask].values
adata.obs["celltype_scanvi"]=adata.obs["celltype_scanvi"].astype("str").astype("category")
labels_key = "celltype_scanvi"

scvi.model.SCVI.setup_anndata(adata, batch_key = batch_key, layer="counts")

vae = scvi.model.SCVI(adata, n_layers=2, n_latent = 30, gene_likelihood="zinb")
vae.train(use_gpu = True)
adata.obsm["X_scVI"]=vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_scVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)


np.unique(adata.obs["celltype_scanvi"], return_counts = True)
lvae = scvi.model.SCANVI.from_scvi_model(
vae,
adata=adata,
unlabeled_category="Unknown",
labels_key="celltype_scanvi",
)
lvae.train( use_gpu=True)
adata.obs["C_scANVI"]=lvae.predict(adata)
adata.obsm["X_scANVI"]=lvae.get_latent_representation(adata)
sc.pp.neighbors(adata, use_rep = "X_scANVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

sc.pl.embedding(adata, basis="X_umap", color=["author_cell_type","C_scANVI"],ncols=1,frameon=False,save=f'_{queryname}_{celltype}_query_ref_scANVI_label_transfer_full.png', palette="tab20")
sc.pl.embedding(adata, basis="X_umap", color=["batch","sample_type","C_scANVI"],ncols=1,frameon=False,save=f'_{queryname}_{celltype}_query_ref_scANVI_label_transfer_sample_type_C_scANVI_full.png', palette="tab20")

ax = sc.pl.umap(
    adata,
    frameon=False,
    show=False,
)
sc.pl.umap(
    adata[adata_ref.n_obs:],
    color=["C_scANVI"],
    frameon=False,
    title="Query predictions",
    ax=ax,
    alpha=0.7,
    save=f'_{queryname}_{celltype}_query_X_scANVI_label_transfer_pred_full.png',
)

ax = sc.pl.umap(
    adata,
    frameon=False,
    show=False,
)
sc.pl.umap(
    adata[adata_ref.n_obs:],
    color=["author_cell_type"],
    frameon=False,
    title="Query observed cell types",
    ax=ax,
    alpha=0.7,
    save=f'_{queryname}_{celltype}_ref_query_X_scANVI_label_transfer_obs_full.png',
)

markers={
"1A":["Pcdh17"],
"1B":["Pcdh10"],
"2":["Ebf1"],
"3A":["Erbb4"],
"3B":["Lrrtm3"],
"4":["Col11a1"],
"5A":["Ptprt"],
"5B":["Chrm2"],
"5C":["Slitrk5"],
"5D":["Kirrel3"],
"6":["Reln"],
"7":["Sox5"],
"8":["Plcb1"],
"9":["Cpne9"],
"RBC":["Prkca"],
}


sc.pl.stacked_violin(

        adata,
        markers,
        groupby='C_scANVI',
        dendrogram=True,
        cmap="Reds",
        swap_axes=True,
        use_raw=True,
        standard_scale="var",
        save=f'_{queryname}_{celltype}_vlnplot.png',
)

sc.pl.dotplot(adata, markers, 'C_scANVI', save=f'_{queryname}_{celltype}_dotplot.png', )


df = adata.obs.groupby(["author_cell_type", "C_scANVI"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]

fig=plt.figure(figsize=(8, 8))
fig=plt.pcolor(conf_mat)
fig= plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
fig= plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
fig=plt.xlabel("Predicted")
fig=plt.ylabel("Observed")
plt.savefig(f'{queryname}_{celltype}_query_ref_scANVI_label_transfer_conf_mat_full.png',format="png",facecolor="w")
conf_mat.to_csv(f'{queryname}_{celltype}_query_ref_scANVI_label_transfer_conf_mat_table_full.csv',header=True, index=True)
adata.write(f'{queryname}_{celltype}_query_ref_label_transfer_full.h5ad')
adata.obs['sample_type']=adata.obs['sample_type'].astype("str").astype("category")

adata_subset = adata[adata.obs['sample_type'].isin(['query'])]
sc.pl.embedding(adata_subset, basis="X_umap", color=["batch","C_scANVI"],ncols=1,frameon=False,save=f'_{queryname}_{celltype}.png', palette="tab20")
sc.pl.dotplot(adata_subset, markers, 'C_scANVI', save=f'_{queryname}_only_{celltype}_dotplot.png', )
adata_subset.write(f'{queryname}_{celltype}.h5ad')
adata_subset.obs['barcode']= adata_subset.obs.index
adata_subset.obs.to_csv(f'{queryname}_{celltype}_obs.txt.gz', sep='\t', index=False)
