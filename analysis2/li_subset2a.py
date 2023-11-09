"""

refine CD8 T cell subsetting 

"""

## import statements
import scanpy as sc
import scanpy.external as sce

## setup
random1 = 1

## file paths
li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
#li_just_T_filename = 'li_T_viz_ready.h5ad'
li_T_wang_genes_filename = 'li_dd_T_wang_genes.h5ad'
li_T_refined_filename = 'li_T_refined.h5ad'

## load the data
adata_li_T = sc.read_h5ad(li_path + li_T_wang_genes_filename)

## process the data

## subset the anndata object -- just the T cells 
#clusters_T = ['3', '4', '5', '6', '7', '8'] ################################
clusters_T = ['0', '2', '3', '4', '5', '6', '7', '8']

adata_li_T_refine = adata_li_T[[x in clusters_T for x in adata_li_T.obs['leiden']]]

## run PCA
print('running PCA...')
sc.tl.pca(adata_li_T_refine, random_state=random1)

## run Harmony batch correction by sample 
print('running Harmony...')
sce.pp.harmony_integrate(adata_li_T_refine, 'patient')

adata_li_T_refine.obsm['X_pca_original'] = adata_li_T_refine.obsm['X_pca']
adata_li_T_refine.obsm['X_pca'] = adata_li_T_refine.obsm['X_pca_harmony']

## compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_li_T_refine, random_state=random1)

## perform Leiden clustering 
print('running Leiden algorithm...')
sc.tl.leiden(adata_li_T_refine, random_state=random1)

## run PAGA trajectory inference 
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_li_T_refine)
sc.pl.paga(adata_li_T_refine, plot=False)

## run UMAP dimensionality reduction 
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_li_T_refine, init_pos='paga', random_state=random1)

## generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

## save to output .h5ad file 
print('writing to output file...')
adata_li_T.write_h5ad(li_path + li_T_refined_filename)
