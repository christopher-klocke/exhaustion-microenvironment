"""
subset to just CD8 T cells

save-as of 'yost_bcc_subset.py'

data from 'remove_rerun_yost_1b.py
"""

## import statements 
import scanpy as sc
import scanpy.external as sce

## setup 
random1 = 1
sample_id = 'sample'

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # file path for local
yost_bcc_just_T_filename = 'data/yost_data/yost_bcc_T.h5ad'
yost_bcc_T_refined_filename = 'data/yost_data/yost_bcc_T_refined.h5ad'

## load the data
adata_yost = sc.read_h5ad(dir + yost_bcc_just_T_filename)

## subset the data
cd8_clusters = ['0', '1', '2', '4', '6', '7', '8', '9', '11']
adata_yost_T = adata_yost[adata_yost.obs['leiden'].isin(cd8_clusters)]

## re-calculate highly variable genes 
del adata_yost_T.var['highly_variable']
adata_yost_T.uns['log1p'] = {'base': None}

sc.pp.highly_variable_genes(adata_yost_T)
adata_yost_T.raw = adata_yost_T

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
adata_yost_T.write_h5ad(dir + yost_bcc_T_refined_filename)
