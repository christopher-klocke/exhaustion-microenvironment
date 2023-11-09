"""

subset to just T cells

"""

## import statements 
import scanpy as sc
import scanpy.external as sce

## setup 
random1 = 1

## file paths 
wang_dir = '/Users/klockec/Documents/data/wang_data/'
wang_viz_ready_filename = 'wang_dd_viz_ready.h5ad'
wang_just_T_filename = 'wang_dd_T_viz_ready.h5ad'

## load the data
adata_wang = sc.read_h5ad(wang_dir + wang_viz_ready_filename)

## process the data

## subset the anndata object -- just the T cells 
clusters_T = ['0', '2', '3', '5', '6'] ## UPDATED 
adata_wang_T = adata_wang[[x in clusters_T for x in adata_wang.obs['leiden']]]

## re-calculate highly variable genes 

"""
troubleshoot -- known error? 

https://github.com/scverse/scanpy/issues/2181

manually set adata.uns['log1p'] to {'base': None}

data are logarithmized though -- check this

"""

del adata_wang_T.var['highly_variable']
adata_wang_T.uns['log1p'] = {'base': None}

sc.pp.highly_variable_genes(adata_wang_T)
adata_wang_T.raw = adata_wang_T ## 

## run PCA
print('running PCA...')
sc.tl.pca(adata_wang_T, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_wang_T, 'batch')

adata_wang_T.obsm['X_pca_original'] = adata_wang_T.obsm['X_pca']
adata_wang_T.obsm['X_pca'] = adata_wang_T.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_wang_T, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
sc.tl.leiden(adata_wang_T, random_state=random1)

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_wang_T)
sc.pl.paga(adata_wang_T, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_wang_T, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_wang_T.write_h5ad(wang_dir + wang_just_T_filename)
