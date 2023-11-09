"""run with 'sc2' conda env """

# import statements
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# setup
adata_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_t_monocle_run4_imputed.h5ad'  # from 'yost_t_impute_updated.py'
umap_embeddings_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/umap_embeddings.csv"  # from 'analysis2/monocle3_r/yost_pt_update.R'
partitions_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/partitions.csv"  # from 'analysis2/monocle3_r/yost_pt_update.R'
#plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/yost_t/updated/'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/yost_t/updated2/'
os.chdir(plots_store)
gene_list = ['TCF7', 'TOX', 'CD8A', 'CD8B', 'PDCD1', 'LAG3', 'TIGIT']
genes1 = ['TCF7', 'TOX', 'CD8B', 'LAG3', 'TIGIT', 'PDCD1']

# load the data
adata = sc.read_h5ad(adata_filename)
df_umap = pd.read_csv(umap_embeddings_filename, names=['CellID', 'UMAP1', 'UMAP2'])
df_umap = df_umap.iloc[1:]
df_umap.index = df_umap['CellID']
df_umap['UMAP1'] = df_umap['UMAP1'].astype(float)
df_umap['UMAP2'] = df_umap['UMAP2'].astype(float)
umap_add = df_umap.to_numpy()
df_partitions = pd.read_csv(partitions_filename)

# add partitions and embeddings
index_check = ['-'.join(x.split('-')[:-1]) for x in adata.obs.index]
df_umap = df_umap.reindex(index_check)
if sum(index_check == df_umap['CellID']) != len(index_check):
    print('ERROR')
    exit()
umap_add = df_umap[['UMAP1', 'UMAP2']].to_numpy()
adata.obsm['X_umap'] = umap_add  # replace UMAP embeddings with those imported from Monocle3

# create plots
#sc.set_figure_params(figsize=(12, 10), fontsize=50)
sc.set_figure_params(figsize=(12, 10), fontsize=50, vector_friendly=False)
size = 20
for gene in genes1:
    sc.pl.umap(adata, color=gene, size=size, show=False, save='_yost_t_' + gene + '.pdf')
sc.pl.umap(adata, color='monocle3_pseudotime', size=size, show=False, save='_yost_t_pt.pdf')
