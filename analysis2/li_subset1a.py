"""

subset to just T cells 

pulling in highly variable genes from wang T cells -- using to define trajectory more clearly 

"""

## import statements 
import scanpy as sc
import scanpy.external as sce

## setup 
random1 = 1

## file paths 
li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
li_viz_ready_filename = 'li_dd_viz_ready.h5ad'
li_just_T_filename = 'li_dd_T_viz_ready.h5ad'
li_T_wang_genes_filename = 'li_dd_T_wang_genes.h5ad'
wang_dir = '/Users/klockec/Documents/data/wang_data/'
wang_just_T_filename = 'wang_dd_T_viz_ready.h5ad'

## load the data
adata_li = sc.read_h5ad(li_path + li_viz_ready_filename)
adata_wang = sc.read_h5ad(wang_dir + wang_just_T_filename)

## process the data

## subset the anndata object -- just the T cells 
#clusters_T = ['3', '4', '5', '6', '7', '8'] ## OLD
clusters_T = ['0', '2', '4', '6', '9', '10'] ## UPDATED 
adata_li_T = adata_li[[x in clusters_T for x in adata_li.obs['leiden']]]

## re-calculate highly variable genes 

"""
troubleshoot -- known error? 

https://github.com/scverse/scanpy/issues/2181

manually set adata.uns['log1p'] to {'base': None}

data are logarithmized though -- check this

"""

del adata_li_T.var['highly_variable']
adata_li_T.uns['log1p'] = {'base': None}

sc.pp.highly_variable_genes(adata_li_T)
adata_li_T.raw = adata_li_T ## 

##############

adata_li_T.var['internal_highly_variable'] = adata_li_T.var['highly_variable'] ## save highly variable genes calculated from Li dataset itself (before replacing with Wang highly variable genes)

del adata_li_T.var['highly_variable']
adata_li_T.uns['log1p'] = {'base': None} ####### CHECK ##############

wang_list = list(adata_wang[:, adata_wang.var['highly_variable']].var.index)

adata_li_T.var['highly_variable'] = [x in wang_list for x in adata_li_T.var.index] ## USE HIGHLY VARIABLE GENES FROM WANG DATASET -- will be used by PCA, fed into downstream neighbor graph construction 

##############

## run PCA
print('running PCA...')
sc.tl.pca(adata_li_T, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_li_T, 'patient')

adata_li_T.obsm['X_pca_original'] = adata_li_T.obsm['X_pca']
adata_li_T.obsm['X_pca'] = adata_li_T.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_li_T, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
sc.tl.leiden(adata_li_T, random_state=random1, resolution=0.95) ####### CHANGING CLUSTERING RESOLUTION

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_li_T)
sc.pl.paga(adata_li_T, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_li_T, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_li_T.write_h5ad(li_path + li_T_wang_genes_filename)
