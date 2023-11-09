"""
generate and save plots for yost dataset

save-as of yost_plots2.py
"""

## import statements
import scanpy as sc
import os

## setup 
random1 = 1
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/yost_full/'
os.chdir(plots_store)

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
yost_main_filename = 'data/yost_data/yost_dd_viz_ready.h5ad'
yost_imputed_filename = 'data/yost_data/yost_dd_imputed.h5ad'
#yost_imputed_add_filename = 'data/yost_data/yost_dd_imputed_add.h5ad'
yost_imputed_add_filename = 'data/yost_data/yost_dd_imputed_add2.h5ad'

## load the data
adata_yost = sc.read_h5ad(dir + yost_main_filename)
adata_yost_imputed = sc.read_h5ad(dir + yost_imputed_filename)
adata_yost_imputed_add = sc.read_h5ad(dir + yost_imputed_add_filename)

## generate figures
main_genes = ['CD8A', 'CD8B', 'GZMK', 'TOX', 'TCF7', 'CD4', 'FOXP3', 'LAG3', 'TIGIT', 'PDCD1', 'CTLA4', 'CD3E']
group1 = ['SSR4', 'IGKC', 'IGHA1', 'FKBP11', 'CD79A', 'MZB1', 'TNFRSF17', 'CD38']
group2 = ['GZMB', 'SEC61B', 'NR3C1', 'CD74', 'RPS12', 'HLA-DRA', 'HLA-DRB1', 'IRF4', 'FUT4', 'TLR7']
group3_1 = ['LYZ', 'TYROBP', 'LGALS1', 'SAT1', 'LMNA', 'HLA-DRB5', 'S100A4', 'HLA-DRA', 'HLA-DRB1', 'CD74', 'RPS12', 'FTL']
group3_2 = ['CST3', 'SIRPA', 'FTH1', 'CSF1R', 'CD14', 'CCR7', 'FLT3', 'CD80', 'CD86', 'TNF', 'XCR1', 'CLEC9A']
group3_3 = ['ITGAM', 'ITGAX', 'ZBTB46', 'SPI1', 'TET2', 'VCAN']
group4 = ['MS4A1', 'BANK1', 'CD79A', 'CD79B', 'LTB', 'CD37', 'CXCR5', 'CD19']
group5 = ['GNLY', 'NKG7', 'CCL5', 'GZMH', 'CST7', 'CCR5', 'EOMES', 'KLRK1']
group6 = ['CD4', 'FOXP3', 'SAT1', 'IL32', 'BATF', 'IL2RA', 'GATA3', 'CD27', 'ADGRE1']
group7 = ['STMN1', 'HMGN2', 'HMGB1', 'CST7']
group8 = ['RGCC', 'MYADM', 'LMNA']
group9 = ['DNAJB1', 'HSPA1B', 'HSPA1A']
groupT_1 = ['SPOCK2', 'TRBV20-1', 'NR3C1', 'S100A4', 'PABPC1', 'B2M', 'IL7R', 'RPS12', 'LTB', 'HAVCR2', 'CXCR3', 'CXCR4']
groupT_2 = ['NFATC1', 'TOX2', 'IFNG', 'SLAMF6', 'PRDM1', 'SELL', 'TBX21', 'TNFRSF9', 'KLRK1', 'PTPRC', 'KLRB1', 'ZAP70']
groupT_3 = ['LCK', 'PLCG1', 'CXCL13']

nk1 = ['CD56', 'CD16', 'CD3', 'NKG7', 'NKG2A', 'KLRB1', 'CD94', 'GNLY', 'KLRC1', 'GZMB', 'NKG2D', 'NCAM1'] ## some genes not in anndata object
#nk1 = ['NKG7', 'KLRB1', 'GNLY', 'KLRC1', 'GZMB', 'NCAM1', 'KLRD1', 'FGFBP2', 'CX3CR1'] ## reduced
#nk1 = ['CD56', 'CD16', 'CD3', 'NKG7', 'NKG2A', 'KLRB1', 'CD94', 'GNLY', 'GZMB', 'NKG2D', 'NCAM1'] ## some genes not in anndata object
nk2 = ['CD45', 'KLRD1', 'FGFBP2', 'FCG3RA', 'CX3CR1'] ## remove genes
#macs1 = ['CD68', 'CD163', 'CD14', 'CD11B'] ## remove genes not in anndata
macs1 = ['CD68', 'CD163', 'CD14']
#plasma_b = ['SDC1', 'CD20'] ## remove missing genes
plasma = ['SDC1']
#b2 = ['CD19', 'CD27', 'CD38', 'IGD', 'CD79A', 'CD37', 'BLNK', 'MS4A1', 'CD24']
b2 = ['CD19', 'CD27', 'CD38', 'CD79A', 'CD37', 'BLNK', 'MS4A1', 'CD24']  # updated -- remove IGD
#plasma2 = ['CD38', 'CD138', 'CD19', 'CD319', 'CD27', 'IGD']
plasma2 = ['CD38', 'CD19', 'CD319', 'CD27']  # updated -- remove IGD, CD138

#sc.pl.umap(adata_yost, show=False, color='leiden', save=('_yost_full_' + 'leiden.pdf')) ##

#sc.pl.umap(adata_yost_imputed, show=False, color=main_genes, save=('_yost_full_' + 'main_genes' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group1, save=('_yost_full_' + 'group1' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group2, save=('_yost_full_' + 'group2' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group3_1, save=('_yost_full_' + 'group3_1' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group3_2, save=('_yost_full_' + 'group3_2' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group3_3, save=('_yost_full_' + 'group3_3' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group4, save=('_yost_full_' + 'group4' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group5, save=('_yost_full_' + 'group5' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group6, save=('_yost_full_' + 'group6' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group7, save=('_yost_full_' + 'group7' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group8, save=('_yost_full_' + 'group8' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=group9, save=('_yost_full_' + 'group9' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=groupT_1, save=('_yost_full_' + 'groupT_1' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=groupT_2, save=('_yost_full_' + 'groupT_2' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed, show=False, color=groupT_3, save=('_yost_full_' + 'groupT_3' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed_add, show=False, color=macs1, save=('_yost_full_' + 'macs1' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed_add, show=False, color=plasma, save=('_yost_full_' + 'plasma_b' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed_add, show=False, color=nk1, save=('_yost_full_' + 'nk1' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed_add, show=False, color=nk2, save=('_yost_full_' + 'nk2' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost_imputed_add, show=False, color=b2, save=('_yost_full_' + 'b2' + '_impute.pdf')) ##
sc.pl.umap(adata_yost_imputed_add, show=False, color=plasma2, save=('_yost_full_' + 'plasma2' + '_impute.pdf')) ##

#sc.pl.umap(adata_yost_imputed, show=False, color=groupT_3, save=('_yost_full_' + 'groupT_3' + '_impute.pdf')) ##
#sc.pl.umap(adata_yost, show=False, color=nk1, save=('_yost_full_' + 'nk1' + '_no_impute.pdf')) ##
#sc.pl.umap(adata_yost, show=False, color=nk2, save=('_yost_full_' + 'nk2' + '_no_impute.pdf')) ##
#sc.pl.umap(adata_yost, show=False, color=macs1, save=('_yost_full_' + 'macs1' + '_no_impute.pdf')) ##
#sc.pl.umap(adata_yost, show=False, color=plasma, save=('_yost_full_' + 'plasma_b' + '_no_impute.pdf')) ##
