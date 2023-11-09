"""
visualize Li data -- just T cell clusters

find starting point for DPT CD8 T cell exhaustion trajectory 
"""

## import statements 
import scanpy as sc
import os

## setup 
random1 = 1
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/li_t/'
os.chdir(plots_store)

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
#li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
#li_just_T_filename = 'li_T_viz_ready.h5ad'
#li_imputed_filename = 'li_T_imputed.h5ad'
#li_T_wang_genes_filename = 'li_T_wang_genes.h5ad'
li_T_wang_genes_filename = 'li_dd_T_wang_genes.h5ad'
#li_imputed_wang_genes_filename = 'li_T_wang_genes_imputed.h5ad'
li_imputed_wang_genes_filename = 'data/li_data/li_T_wang_genes_imputed.h5ad'

## load the data
#adata_li_T = sc.read_h5ad(li_path + li_T_wang_genes_filename)
adata_li_T = sc.read_h5ad(dir + li_imputed_wang_genes_filename)

## generate plots
#genes_plot = ['CD8A', 'CD8B', 'PDCD1', 'TOX', 'LAG3', 'TIGIT', 'TCF7', 'CD3E', 'CD4', 'FOXP3', 'GNLY', 'NKG7']
gene_list1 = ['CD3E', 'CD3G', 'CD4', 'FOXP3', 'CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'LAG3', 'TIGIT']
gene_list2 = ['PDCD1', 'NKG7', 'GNLY', 'TYROBP', 'KLRK1', 'CST7', 'ITGAM', 'NCAM1']

#sc.pl.umap(adata_li_T, color='leiden', show=False, save='_li_t_leiden_d95.pdf')
#sc.pl.umap(adata_li_T, color='leiden')
#sc.pl.umap(adata_li_T, color=gene_list1, show=False, save='_li_t_genes1.pdf')
#sc.pl.umap(adata_li_T, color=gene_list1)
#sc.pl.umap(adata_li_T, color=gene_list2, show=False, save='_li_t_genes2.pdf')
#sc.pl.umap(adata_li_T, color=gene_list2)
#sc.pl.umap(adata_li_T, color='leiden')

sc.pl.umap(adata_li_T, color='TCF7', show=False, save='_li_t_TCF7.pdf')
sc.pl.umap(adata_li_T, color='TOX', show=False, save='_li_t_TOX.pdf')
sc.pl.umap(adata_li_T, color='LAG3', show=False, save='_li_t_LAG3.pdf')
sc.pl.umap(adata_li_T, color='TIGIT', show=False, save='_li_t_TIGIT.pdf')
sc.pl.umap(adata_li_T, color='PDCD1', show=False, save='_li_t_PDCD1.pdf')
sc.pl.umap(adata_li_T, color='CD8B', show=False, save='_li_t_CD8B.pdf')
