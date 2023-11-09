"""
save-as of 'remove_rerun_li.py'

run with 'sc2' conda env
"""

## import statements
import scanpy as sc
import scanpy.external as sce

## setup
random1 = 1
sample_id = 'sample'

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
yost_path = 'data/yost_data/'  # updated for local
yost_doublets_adata_filename = 'yost_doublets_adata.h5ad'
yost_output_filename = 'yost_dd_viz_ready.h5ad'

## load the data
adata_yost = sc.read_h5ad(dir + yost_path + yost_doublets_adata_filename)

## remove the doublets
adata_yost = adata_yost[adata_yost.obs['doublet'] == 0]

## re-process the data

## compute highly variable genes
del adata_yost.var['highly_variable']
adata_yost.uns['log1p'] = {'base': None}

sc.pp.highly_variable_genes(adata_yost)
adata_yost.raw = adata_yost

## run PCA
print('running PCA...')
sc.tl.pca(adata_yost, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_yost, sample_id)

adata_yost.obsm['X_pca_original'] = adata_yost.obsm['X_pca']
adata_yost.obsm['X_pca'] = adata_yost.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_yost, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
sc.tl.leiden(adata_yost, random_state=random1)

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_yost)
sc.pl.paga(adata_yost, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_yost, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_yost.write_h5ad(dir + yost_path + yost_output_filename)
