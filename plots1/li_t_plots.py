"""
visualize Li data -- just T cell clusters

find starting point for DPT CD8 T cell exhaustion trajectory 
"""

## import statements 
import scanpy as sc
import os

## setup 
random1 = 1
plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_li_t'
os.chdir(plots_store)

## file paths 
li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
li_just_T_filename = 'li_T_viz_ready.h5ad'
li_imputed_filename = 'li_T_imputed.h5ad'

## load the data
adata_li_T = sc.read_h5ad(li_path + li_imputed_filename)

## generate plots

#genes_plot = ['CD8A', 'CD8B', 'PDCD1', 'TOX', 'LAG3', 'TIGIT', 'TCF7', 'CD3E', 'CD4', 'FOXP3', 'GNLY', 'NKG7']
gene_list1 = ['CD3E', 'CD3G', 'CD4', 'FOXP3', 'CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'LAG3', 'TIGIT']
gene_list2 = ['PDCD1', 'NKG7', 'GNLY', 'TYROBP', 'KLRK1', 'CST7', 'ITGAM', 'NCAM1']

sc.pl.umap(adata_li_T, color='leiden', save='_li_t_leiden_d95.pdf')
#sc.pl.umap(adata_li_T, color='leiden')
sc.pl.umap(adata_li_T, color=gene_list1, save='_li_t_genes1.pdf')
#sc.pl.umap(adata_li_T, color=gene_list1)
sc.pl.umap(adata_li_T, color=gene_list2, save='_li_t_genes2.pdf')
#sc.pl.umap(adata_li_T, color=gene_list2)
#sc.pl.umap(adata_li_T, color='leiden')
