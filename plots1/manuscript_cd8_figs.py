"""
imputed files from 'manuscript_cd8_impute.py'
"""
#  import statements
import os
import scanpy as sc

#  setup
dir_name = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
li_imputed_filename = 'data/li_data/li_t_monocle1_imputed.h5ad'
yost_imputed_filename = 'data/yost_data/yost_t_monocle1_v2_imputed.h5ad'
wang_imputed_filename = 'data/wang_data/wang_t_monocle1_imputed.h5ad'
genes1 = ['TCF7', 'TOX', 'LAG3', 'CD8B', 'PDCD1', 'TIGIT']
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/cd8_all/'
pt_name = 'monocle3_pseudotime'

os.chdir(plots_store)

#  load the data
adata_li = sc.read_h5ad(dir_name + li_imputed_filename)
adata_yost = sc.read_h5ad(dir_name + yost_imputed_filename)
adata_wang = sc.read_h5ad(dir_name + wang_imputed_filename)

#  generate visualizations
sc.pl.umap(adata_li, color=pt_name, save='_li_pt.pdf')
sc.pl.umap(adata_yost, color=pt_name, save='_yost_pt.pdf')
sc.pl.umap(adata_wang, color=pt_name, save='_wang_pt.pdf')

for gene in genes1:
    sc.pl.umap(adata_li, color=gene, save='_li_'+gene+'_imputed.pdf')
    sc.pl.umap(adata_yost, color=gene, save='_yost_'+gene+'_imputed.pdf')
    sc.pl.umap(adata_wang, color=gene, save='_wang_'+gene+'_imputed.pdf')

