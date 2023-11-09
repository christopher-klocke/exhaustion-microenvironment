"""
perform imputation with the MAGIC algorithm

runs in ~26min for all 3 datasets

this script is save-as of impute.py

marking as v3 to be in accordance with li_plots3.py and wang_plots3.py
"""

## import statements
import scanpy as sc
import magic

## setup
random1 = 1

## file paths
#li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
#yost_path = '/Users/klockec/Documents/data/yost_data/'
#wang_path = '/Users/klockec/Documents/data/wang_data/'
#li_filename = 'li_dd_viz_ready.h5ad'
#yost_filename = 'yost_dd_viz_ready.h5ad'
#wang_filename = 'wang_dd_viz_ready.h5ad'
#li_imputed_filename = 'li_dd_imputed.h5ad'
#yost_imputed_filename = 'yost_dd_imputed.h5ad'
#wang_imputed_filename = 'wang_dd_imputed.h5ad'

## file paths (updated for local)
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
yost_filename = 'data/yost_data/yost_dd_viz_ready.h5ad'
#yost_imputed_filename = 'data/yost_data/yost_dd_imputed_add.h5ad'
yost_imputed_filename = 'data/yost_data/yost_dd_imputed_add2.h5ad'

## load the data
#adata_li = sc.read_h5ad(li_path + li_filename)
adata_yost = sc.read_h5ad(dir + yost_filename)
#adata_wang = sc.read_h5ad(wang_path + wang_filename)

## list of genes for imputation
#gene_list = ['CD56', 'CD16', 'CD3', 'NKG7', 'NKG2A', 'KLRB1', 'CD94', 'GNLY', 'KLRC1', 'GZMB', 'NKG2D', 'NCAM1', 'CD45', 'KLRD1', 'FGFBP2', 'FCG3RA', 'CX3CR1', 'CD68', 'CD163', 'CD14', 'CD11B', 'SDC1', 'CD20']
gene_list = ['CD19', 'CD27', 'CD38', 'IGD', 'CD79A', 'CD37', 'BLNK', 'MS4A1', 'CD24', 'CD138', 'CD319']

## impute with MAGIC algorithm

"""
### li data
X_li = adata_li
magic_operator = magic.MAGIC()
X_li_magic = magic_operator.fit_transform(X_li, genes=gene_list)
X_li_magic.uns['umap'] = adata_li.uns['umap']
X_li_magic.obsm['X_umap'] = adata_li.obsm['X_umap']
"""

### yost data
X_yost = adata_yost
magic_operator = magic.MAGIC()
X_yost_magic = magic_operator.fit_transform(X_yost, genes=gene_list)
X_yost_magic.uns['umap'] = adata_yost.uns['umap']
X_yost_magic.obsm['X_umap'] = adata_yost.obsm['X_umap']

""""
### wang data
X_wang = adata_wang
#magic_operator = magic.MAGIC()
X_wang_magic = magic_operator.fit_transform(X_wang, genes=gene_list)
X_wang_magic.uns['umap'] = adata_wang.uns['umap']
X_wang_magic.obsm['X_umap'] = adata_wang.obsm['X_umap']
"""

## write to output files
#X_li_magic.write_h5ad(li_path + li_imputed_filename)
X_yost_magic.write_h5ad(dir + yost_imputed_filename)
#X_wang_magic.write_h5ad(wang_path + wang_imputed_filename)
