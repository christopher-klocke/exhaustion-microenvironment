"""
visualize Wang data -- just T cell clusters 

find starting point for DPT CD8 T cell exhaustion trajectory 
"""

## import statements 
import scanpy as sc
import os

## setup 
random1 = 1
#plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_wang_runA'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/wang_t/'
os.chdir(plots_store)

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
#wang_dir = '/Users/klockec/Documents/data/wang_data/'
#wang_just_T_filename = 'wang_T_viz_ready.h5ad'
#wang_just_T_filename = 'wang_dd_T_viz_ready.h5ad'
wang_just_T_filename = 'data/wang_data/wang_dd_T_viz_ready.h5ad'  # updated for local location
wang_just_T_imputed_filename = 'data/wang_data/wang_T_imputed.h5ad'

## load the data
#adata_wang_T = sc.read_h5ad(wang_dir + wang_just_T_filename)
#adata_wang_T = sc.read_h5ad(dir + wang_just_T_filename)
adata_wang_T_imputed = sc.read_h5ad(dir + wang_just_T_imputed_filename)

## generate plots
#genes_plot = ['CD8A', 'CD8B', 'PDCD1', 'TOX', 'LAG3', 'TIGIT', 'TCF7', 'CD3E']
#sc.pl.umap(adata_wang_T, color=genes_plot)
#sc.pl.umap(adata_wang_T, color='leiden')

#gene_list1 = ['CD3E', 'CD3G', 'CD4', 'FOXP3', 'CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'LAG3', 'TIGIT']
#gene_list2 = ['PDCD1', 'NKG7', 'GNLY', 'TYROBP', 'KLRK1', 'CST7', 'ITGAM', 'NCAM1']

#sc.pl.umap(adata_wang_T, color='leiden', show=False, save='_wang_t_leiden.pdf')
#sc.pl.umap(adata_wang_T, color=gene_list1, show=False, save='_wang_t_genes1.pdf')
#sc.pl.umap(adata_wang_T, color=gene_list2, show=False, save='_wang_t_genes2.pdf')

#sc.pl.umap(adata_wang_T_imputed, color='TCF7', show=False, save='_wang_t_TCF7.pdf')
#sc.pl.umap(adata_wang_T_imputed, color='TOX', show=False, save='_wang_t_TOX.pdf')
#sc.pl.umap(adata_wang_T_imputed, color='LAG3', show=False, save='_wang_t_LAG3.pdf')
#sc.pl.umap(adata_wang_T_imputed, color='TIGIT', show=False, save='_wang_t_TIGIT.pdf')
#sc.pl.umap(adata_wang_T_imputed, color='PDCD1', show=False, save='_wang_t_PDCD1.pdf')

sc.pl.umap(adata_wang_T_imputed, color='CD3E', show=False, save='_wang_t_CD3E.pdf')
sc.pl.umap(adata_wang_T_imputed, color='CD8A', show=False, save='_wang_t_CD8A.pdf')
sc.pl.umap(adata_wang_T_imputed, color='CD8B', show=False, save='_wang_t_CD8B.pdf')
