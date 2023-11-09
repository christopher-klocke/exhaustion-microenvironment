"""
subset to just T cells

save-as of 'yost_subset1.py'

subset again
"""

## import statements 
import scanpy as sc
import scanpy.external as sce

## setup 
random1 = 1
sample_id = 'sample'

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
#yost_filename = 'data/yost_data/yost_dd_viz_ready.h5ad'
yost_just_T_filename1 = 'data/yost_data/yost_dd_T_viz_ready_v2.h5ad'
yost_just_T_filename2 = 'data/yost_data/yost_dd_T_viz_ready_v3.h5ad'

## load the data
adata_yost = sc.read_h5ad(dir + yost_just_T_filename1)

## process the data

## subset the anndata object -- just the CD8 T cells
#clusters_T = ['0', '2', '3', '5', '6']
#adata_yost_T = adata_yost[[x in clusters_T for x in adata_yost.obs['leiden']]]
#['CD8_ex', 'CD8_mem', 'Treg', 'Tfh', 'Naive', 'CD8_act', 'CD8_naive', 'CD8_ex_act', 'Th17', 'CD8_eff', 'CD8_ex_T_cells', 'CD8_mem_T_cells', 'CD4_T_cells', 'Tregs', 'Tcell_prolif', 'CD8_act_T_cells', 'Macrophages', 'Plasma_cells', 'DCs', 'pDCs', 'B_cells_1', 'NK_cells', 'B_cells_2']

#clusters_T = ['CD8_ex', 'CD8_mem', 'CD8_act', 'CD8_naive', 'CD8_ex_act', 'CD8_eff', 'CD8_ex_T_cells', 'CD8_mem_T_cells', 'CD8_act_T_cells']
#adata_yost_T = adata_yost[[x in clusters_T for x in adata_yost.obs['cluster_original']]]

clusters_T = ['3', '4', '5', '8', '9', '10', '11']
adata_yost_T = adata_yost[[x in clusters_T for x in adata_yost.obs['leiden']]]

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
adata_yost_T.write_h5ad(dir + yost_just_T_filename2)
