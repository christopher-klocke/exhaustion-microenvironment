"""
save-as of 'remove_rerun_yost.py'

run with 'sc2' conda env

input files from 'doublet_yost_dd2b.py'
"""

## import statements
import scanpy as sc
import scanpy.external as sce

## setup
random1 = 1
sample_id = 'sample'

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # file path for local
yost_bcc_doublets_adata_filename = 'data/yost_data/yost_bcc_doublets_adata.h5ad'
yost_scc_doublets_adata_filename = 'data/yost_data/yost_scc_doublets_adata.h5ad'
yost_bcc_full_filename = 'data/yost_data/yost_bcc_full_processed_dd.h5ad'
yost_scc_full_filename = 'data/yost_data/yost_scc_full_processed_dd.h5ad'

## load the data
adata_yost_bcc = sc.read_h5ad(dir + yost_bcc_doublets_adata_filename)
adata_yost_scc = sc.read_h5ad(dir + yost_scc_doublets_adata_filename)

## remove the doublets
adata_yost_bcc = adata_yost_bcc[adata_yost_bcc.obs['doublet'] == 0]
adata_yost_scc = adata_yost_scc[adata_yost_scc.obs['doublet'] == 0]

## re-process the data

## compute highly variable genes
del adata_yost_bcc.var['highly_variable']
adata_yost_bcc.uns['log1p'] = {'base': None}
sc.pp.highly_variable_genes(adata_yost_bcc)
adata_yost_bcc.raw = adata_yost_bcc

del adata_yost_scc.var['highly_variable']
adata_yost_scc.uns['log1p'] = {'base': None}
sc.pp.highly_variable_genes(adata_yost_scc)
adata_yost_scc.raw = adata_yost_scc

## run PCA
print('running PCA...')
sc.tl.pca(adata_yost_bcc, random_state=random1)
sc.tl.pca(adata_yost_scc, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_yost_bcc, sample_id)
adata_yost_bcc.obsm['X_pca_original'] = adata_yost_bcc.obsm['X_pca']
adata_yost_bcc.obsm['X_pca'] = adata_yost_bcc.obsm['X_pca_harmony']

sce.pp.harmony_integrate(adata_yost_scc, sample_id)
adata_yost_scc.obsm['X_pca_original'] = adata_yost_scc.obsm['X_pca']
adata_yost_scc.obsm['X_pca'] = adata_yost_scc.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_yost_bcc, random_state=random1)
sc.pp.neighbors(adata_yost_scc, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
sc.tl.leiden(adata_yost_bcc, random_state=random1)
sc.tl.leiden(adata_yost_scc, random_state=random1)

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_yost_bcc)
sc.pl.paga(adata_yost_bcc, plot=False)
sc.tl.paga(adata_yost_scc)
sc.pl.paga(adata_yost_scc, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_yost_bcc, init_pos='paga', random_state=random1)
sc.tl.umap(adata_yost_scc, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_yost_bcc.write_h5ad(dir + yost_bcc_full_filename)
adata_yost_scc.write_h5ad(dir + yost_scc_full_filename)

