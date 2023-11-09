"""
subset to just T cells

pulling in highly variable genes from wang T cells -- using to define trajectory more clearly

save-as of li_subset1a.py
"""

## import statements 
import scanpy as sc
import scanpy.external as sce

## setup 
random1 = 1

## setup
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
yost_T_filename = 'data/yost_data/yost_dd_T_viz_ready_v4.h5ad'
yost_T_wang_genes_filename = 'data/yost_data/yost_T_wang_genes_v1.h5ad'
wang_T_filename = 'data/wang_data/wang_dd_T_viz_ready.h5ad'
sample_id = 'sample'

## load the data
adata_yost_T = sc.read_h5ad(dir + yost_T_filename)
adata_wang_T = sc.read_h5ad(dir + wang_T_filename)

## process the data
adata_yost_T.var['internal_highly_variable'] = adata_yost_T.var['highly_variable'] ## save highly variable genes calculated from dataset itself (before replacing with Wang highly variable genes)
del adata_yost_T.var['highly_variable']
adata_yost_T.uns['log1p'] = {'base': None}
wang_list = list(adata_wang_T[:, adata_wang_T.var['highly_variable']].var.index)
adata_yost_T.var['highly_variable'] = [x in wang_list for x in adata_yost_T.var.index] ## USE HIGHLY VARIABLE GENES FROM WANG DATASET -- will be used by PCA, fed into downstream neighbor graph construction

## run PCA
print('running PCA...')
sc.tl.pca(adata_yost_T, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_yost_T, sample_id)

adata_yost_T.obsm['X_pca_original'] = adata_yost_T.obsm['X_pca']
adata_yost_T.obsm['X_pca'] = adata_yost_T.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_yost_T, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
#sc.tl.leiden(adata_yost_T, random_state=random1, resolution=0.95) ####### CHANGING CLUSTERING RESOLUTION
sc.tl.leiden(adata_yost_T, random_state=random1)

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_yost_T)
sc.pl.paga(adata_yost_T, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_yost_T, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_yost_T.write_h5ad(dir + yost_T_wang_genes_filename)
