"""
visualize Yost data -- just T cell clusters

find starting point for monocle3 CD8 T cell exhaustion trajectory

save-as of 'li_t_plots.py'
"""

## import statements 

import scanpy as sc
import os

## setup 
random1 = 1
#plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/img_yost_t/'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/yost_t/'
os.chdir(plots_store)

## file paths 
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
#yost_just_T_filename = 'data/yost_data/yost_dd_T_viz_ready.h5ad'
#yost_imputed_filename = 'data/yost_data/yost_T_imputed.h5ad'
#yost_imputed_filename = 'data/yost_data/yost_T_imputed_v2.h5ad'
#yost_imputed_filename = 'data/yost_data/yost_T_imputed_v4.h5ad'
yost_imputed_filename = 'data/yost_data/yost_T_imputed_wang_v1.h5ad'

## load the data
adata_yost_T = sc.read_h5ad(dir + yost_imputed_filename)

## generate plots

#genes_plot = ['CD8A', 'CD8B', 'PDCD1', 'TOX', 'LAG3', 'TIGIT', 'TCF7', 'CD3E', 'CD4', 'FOXP3', 'GNLY', 'NKG7']
gene_list1 = ['CD3E', 'CD3G', 'CD4', 'FOXP3', 'CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'LAG3', 'TIGIT']
gene_list2 = ['PDCD1', 'NKG7', 'GNLY', 'TYROBP', 'KLRK1', 'CST7', 'ITGAM', 'NCAM1']

sc.pl.umap(adata_yost_T, color='leiden', save='_yost_t_leiden_wang_v1.pdf')
#sc.pl.umap(adata_yost_T, color='leiden')

sc.pl.umap(adata_yost_T, color=gene_list1, save='_yost_t_genes1_wang_v1.pdf')
#sc.pl.umap(adata_yost_T, color=gene_list1)

sc.pl.umap(adata_yost_T, color=gene_list2, save='_yost_t_genes2_wang_v1.pdf')
#sc.pl.umap(adata_yost_T, color=gene_list2)

sc.pl.umap(adata_yost_T, color='cluster_original', save='_yost_t_cluster_original_wang_v1.pdf')
#sc.pl.umap(adata_yost_T, color='cluster_original')
