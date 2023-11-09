"""

generate and save plots for li dataset

second pass

"""

## import statements

import scanpy as sc
import os

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img2'

os.chdir(plots_store)

## file paths

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

li_main_filename = 'li_viz_ready.h5ad'

li_imputed_filename = 'li_imputed.h5ad'

## load the data

adata_li = sc.read_h5ad(li_path + li_main_filename)

adata_li_imputed = sc.read_h5ad(li_path + li_imputed_filename)

######################

cluster10 = ['IGHG3', 'IGHG1', 'SSR4', 'MZB1', 'CD79A', 'IGJ', 'IRF4', 'PRDM1', 'CD27', 'TNFRSF17']

cluster2 =  ['CD22', 'HLA-DRA', 'CD74', 'MS4A1', 'CD37', 'BANK1', 'LTB', 'CD79A', 'CD79B', 'CXCR5', 'CD19']

cluster12 = ['GZMB', 'CD4', 'TCF4', 'SELL', 'TLR7', 'TLR9', 'HLA-DRA', 'CXCR3', 'IRF4', 'TYROBP', 'IGJ']

cluster9 = ['VIM', 'CD63', 'MIF', 'FTH1']

cluster1_11_p1 = ['VIM', 'IFI30', 'TYROBP', 'FTL', 'CST3', 'CSF1R', 'FTH1', 'SIRPA', 'PTPRC', 'SPI1']

cluster1_11_p2 = ['NFATC1', 'HLA-DRA', 'LYZ', 'CD74', 'CD14', 'CD86', 'VCAN', 'ITGAM', 'ITGAX', 'TET2']

t_p1 = ['CD3E', 'CD8A', 'CD8B', 'CD4', 'FOXP3', 'GZMK', 'GZMB', 'TOX' , 'TCF7', 'LAG3', 'TIGIT', 'PDCD1']

t_p2 = ['TRAC', 'GZMA', 'SYNE2', 'CD2', 'TRBC2', 'CD3G', 'GNLY', 'CCL5', 'NKG7', 'IL32', 'TYROBP', 'IL7R']

t_p3 = ['CTLA4', 'RPS12', 'LTB', 'GZMH', 'CST7', 'IL2RA', 'BATF', 'CX3CR1', 'CXCR3', 'HAVCR2', 'TNFRSF9', 'TBX21']

t_p4 = ['EOMES', 'PRDM1', 'SLAMF6', 'CD101', 'CCR5', 'IFNG', 'TOX2', 'NCAM1', 'PTPRC', 'KLRK1', 'GATA3', 'TNF']

t_p5 = ['KLRB1', 'CD27', 'CD38', 'ZAP70', 'ITGAM', 'PLCG1', 'LCK']

sc.pl.umap(adata_li_imputed, show=False, color=cluster10, save=('_li_full_' + 'cluster10' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=cluster2, save=('_li_full_' + 'cluster2' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=cluster12, save=('_li_full_' + 'cluster12' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=cluster9, save=('_li_full_' + 'cluster9' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=cluster1_11_p1, save=('_li_full_' + 'cluster1_11_p1' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=cluster1_11_p2, save=('_li_full_' + 'cluster1_11_p2' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=t_p1, save=('_li_full_' + 't_p1' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=t_p2, save=('_li_full_' + 't_p2' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=t_p3, save=('_li_full_' + 't_p3' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=t_p4, save=('_li_full_' + 't_p4' + '_impute.pdf'))

sc.pl.umap(adata_li_imputed, show=False, color=t_p5, save=('_li_full_' + 't_p5' + '_impute.pdf'))

