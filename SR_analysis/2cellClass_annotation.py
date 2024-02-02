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

query=sys.argv[1]
celltype=sys.argv[2]
label=sys.argv[3]
queryname=sys.argv[4]

scvi.settings.seed = 0
dataset1=sc.read_10x_mtx(query)
dataset1.obs["sample"]=queryname
dataset1.obs.index=queryname+"-"+dataset1.obs.index

if celltype == "major":
	adata_query=dataset1
else:
	adata_query=dataset1[(dataset1.obs[label] == celltype )]

del adata_query.raw

batch_key = "sample"
labels_key = ""
adata_query.obs[batch_key] = adata_query.obs[batch_key].astype("str").astype("category")
adata_query.layers["counts"] = adata_query.X.copy()
sc.pp.normalize_total(adata_query,target_sum=1e4)
sc.pp.log1p(adata_query)
adata_query.raw = adata_query
try: 
	sc.pp.highly_variable_genes(
		adata_query,
		n_top_genes=2000,
		batch_key = batch_key,
		subset=True,
	)
except:
	print("An exception occurred! Fix the batch or run without batch now.")
	sc.pp.highly_variable_genes(
		adata_query, n_top_genes=2000, subset=True
	)

scvi.settings.seed = 0 
if labels_key == "":
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, layer="counts")

else:
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, labels_key = labels_key)

vae = scvi.model.SCVI(adata_query, n_layers=2, n_latent=30, gene_likelihood="zinb")
vae.train(use_gpu=True)
adata_query.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata_query, use_rep = "X_scVI")
sc.tl.leiden(adata_query, resolution=0.3)
sc.tl.umap(adata_query)
print(adata_query.obs.leiden) 
author_cell_type = adata_query.obs.leiden.astype(str).astype("category")
adata_query.obs["author_cell_type"] = author_cell_type
sc.pl.embedding(adata_query, basis="X_umap", color=["author_cell_type"],ncols=1,frameon=False,save=f'_{queryname}_{celltype}_query_after_label_transfer_reCluster.png', palette="tab20")

adata_query.write(f'{queryname}_{celltype}.h5ad')

markers={"AC":["Gad2","Pax6","Slc6a1","Slc6a9","Slc32a1","Tfap2a"], "BC":["Crx","Grik1","Grm6","Isl1","Otx2","Vsx2"],"Cone":["Adrb1","Casp7","Cngb3"],"HC":["Adra2a","Lima1","Onecut1","Onecut2"],"MG":["Id3","Kdr","Rdh10","Vim"],"RBC":["Prkca","Strip2","Tpbg","Vstm2b"],"RGC":["Pou4f1","Pou4f2","Rbpms","Slc17a6"],"Rod":["Cngb1","Nt5e","Pde6a","Reep6"]}

sc.pl.stacked_violin(

        adata_query,
        markers,
        groupby='author_cell_type',
        dendrogram=True,
        cmap="Reds",
        swap_axes=True,
        use_raw=True,
        standard_scale="var",
        save=f'_{queryname}_{celltype}_vlnplot.png',

)

#From the dotplot, we can know which cluster is which cell type.
sc.pl.dotplot(adata_query, markers, 'author_cell_type', dendrogram=True, save=f'_{queryname}_{celltype}_dotplot.png', )

#assign the cell type name to clusters
# create a dictionary to map cluster to annotation label
#cluster2annotation = {
#     '0': 'Rod', 
#     '1': 'Cone', 
#     '2': 'MG',
#     '3': 'BC', 
#     '4': 'BC',  
#     '5': 'AC',
#     '6': 'BC', 
#     '7': 'BC', 
#     '8': 'BC', 
#     '9': 'BC', 
#     '10': 'BC', 
#     '11': 'MG'
#} 

adata_query.obs['majortype'] = adata_query.obs['author_cell_type'].map(cluster2annotation).astype('category')
 
adata_query.obs['barcode']= adata_query.obs.index
adata_query.obs.to_csv(f'{queryname}_{celltype}_obs.txt.gz', sep='\t', index=False)