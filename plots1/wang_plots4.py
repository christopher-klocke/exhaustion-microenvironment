"""
generate and save plots for wang dataset

save-as of wang_plots3.py
"""

## import statements

import scanpy as sc
import os

## file paths

#dir = '/Users/klockec/OneDrive - Oregon Health & Science University/'
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
wang_path = 'data/wang_data/'  # updated for local
wang_main_filename = 'wang_dd_viz_ready.h5ad'
wang_imputed_filename = 'wang_dd_imputed.h5ad'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/part1/'

## setup 

random1 = 1
os.chdir(plots_store)

## load the data

adata_wang = sc.read_h5ad(dir + wang_path + wang_main_filename)
adata_wang_imputed = sc.read_h5ad(dir + wang_path + wang_imputed_filename)

## plot figures

top_genes = ['CD3E', 'CD8A', 'CD8B', 'TCF7', 'TOX', 'LAG3', 'TIGIT', 'CD4']
sc.pl.umap(adata_wang_imputed, show=False, color=top_genes, save=('_wang_full_' + 'top_genes' + '_impute.pdf'))
#sc.pl.umap(adata_wang_imputed, show=False, color=top_genes)
sc.pl.umap(adata_wang, show=False, color='leiden', save=('_wang_full_leiden.pdf'))
#sc.pl.umap(adata_wang, show=False, color='leiden')

gene_list3 = ['IGJ', 'LTB', 'LYZ']
sc.pl.umap(adata_wang, color=gene_list3, show=False, save='_wang_full_genes3.pdf')

clusters0_5_11_p1 = ['CD3E', 'CD8A', 'CD8B', 'CD4', 'TCF7', 'FOXP3', 'TNFAIP3', 'IL7R', 'RPS12', 'LTB', 'CCR7', 'SELL']
clusters0_5_11_p2 = ['NFATC1', 'GATA3', 'PLCG1']
clusters2_3_6_p1 = ['TOX', 'LAG3', 'TIGIT', 'PDCD1', 'GZMB', 'GZMK', 'FCER1G', 'GNLY', 'CCL5', 'NKG7', 'B2M', 'TYROBP']
clusters2_3_6_p2 = ['IL32', 'GZMH', 'CST7', 'BATF', 'CX3CR1', 'CCR5', 'IFNG', 'SLAMF6', 'TOX2', 'PRDM1', 'EOMES', 'TBX21']
clusters2_3_6_p3 = ['TNFRSF9', 'NCAM1', 'KLRK1', 'PTPRC', 'KLRB1', 'TNF', 'LCK', 'ZAP70']
cluster7_p1 = ['HSP90B1', 'FCER1G', 'SSR4', 'CD79A', 'MZB1', 'IGJ', 'IRF4', 'TOX2', 'PRDM1', 'FLT3', 'TNFRSF17', 'CD27']
cluster7_p2 = ['TLR7', 'CD38']
clusters4_9_p1 = ['CD4', 'LST1', 'AIF1', 'COTL1', 'CTSS', 'HLA-DPA1', 'HLA-DPB1', 'FCER1G', 'TYROBP', 'CD74', 'HLA-DRA', 'LYZ']
clusters4_9_p2 = ['FTL', 'HLA-DRB1', 'CST3', 'FTH1', 'CD14', 'CSF1R', 'CX3CR1', 'HAVCR2', 'CD101', 'SELL', 'NFATC1', 'SIRPA']
clusters4_9_p3 = ['FLT3', 'VCAN', 'FUT4', 'CD86', 'CLEC9A', 'MERTK', 'ITGAM', 'ITGAX', 'ZBTB46', 'CEBPA', 'SPI1', 'TET2']
clusters1_8_p1 = ['HLA-DPA1', 'HLA-DPB1', 'TOX2', 'SLAMF6', 'TBX21', 'CCR7', 'SELL', 'NFATC1', 'CD86', 'CD80', 'CD19', 'TET2']
clusters1_8_p2 = ['CD74', 'HLA-DRA', 'MS4A1', 'CD79A', 'CD79B', 'HLA-DRB1', 'CD37', 'BANK1', 'LTB', 'CXCR4', 'CXCR5']
minor1 = ['PEAR1', 'SDPR', 'HIST1H2AC', 'NAP1L1']
minor2 = ['ALAS2', 'HBB', 'HBA2']

sc.pl.umap(adata_wang_imputed, show=False, color=clusters0_5_11_p1, save=('_wang_full_' + 'clusters0_5_11_p1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters0_5_11_p2, save=('_wang_full_' + 'clusters0_5_11_p2' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters2_3_6_p1, save=('_wang_full_' + 'clusters2_3_6_p1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters2_3_6_p2, save=('_wang_full_' + 'clusters2_3_6_p2' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters2_3_6_p3, save=('_wang_full_' + 'clusters2_3_6_p3' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=cluster7_p1, save=('_wang_full_' + 'cluster7_p1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=cluster7_p2, save=('_wang_full_' + 'cluster7_p2' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters4_9_p1, save=('_wang_full_' + 'clusters4_9_p1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters4_9_p2, save=('_wang_full_' + 'clusters4_9_p2' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters4_9_p3, save=('_wang_full_' + 'clusters4_9_p3' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters1_8_p1, save=('_wang_full_' + 'clusters1_8_p1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=clusters1_8_p2, save=('_wang_full_' + 'clusters1_8_p2' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=minor1, save=('_wang_full_' + 'minor1' + '_impute.pdf'))
sc.pl.umap(adata_wang_imputed, show=False, color=minor2, save=('_wang_full_' + 'minor2' + '_impute.pdf'))
