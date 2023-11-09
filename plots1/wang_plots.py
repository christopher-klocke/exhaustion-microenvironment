"""

generate and save plots for wang dataset

first, generate many plots to look through genes and see which ones are of interest

then, run script again and generate a smaller number of select figures that will be exported

"""

## import statements

import scanpy as sc
import os

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img'

os.chdir(plots_store)

## file paths

wang_path = '/Users/klockec/Documents/data/wang_data/'

wang_main_filename = 'wang_viz_ready.h5ad'

wang_imputed_filename = 'wang_imputed.h5ad'

## load the data

adata_wang = sc.read_h5ad(wang_path + wang_main_filename)

adata_wang_imputed = sc.read_h5ad(wang_path + wang_imputed_filename)

######################

## first, generate plots by cluster and by patient sample

sc.pl.umap(adata_wang, show=False, color='leiden', save=('_pass1_' + 'wang_full_leiden.pdf'))

sc.pl.umap(adata_wang, show=False, color='batch', save=('_pass1_' + 'wang_full_patient.pdf'))

## next, plot subsets of gene list (too many genes to fit well in one plot) -- with and without imputation 

main_genes = ['CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'CD4', 'FOXP3', 'LAG3', 'TIGIT', 'PDCD1', 'CTLA4', 'CD3E']

dataset_specific1 = ['CTSS', 'HLA-DPB1', 'HLA-DPA1', 'TNFAIP3', 'HSP90B1', 'LST1', 'AIF1', 'COTL1']

dataset_specific2 = ['FCER1G', 'HBA2', 'HBB', 'ALAS2', 'SDPR', 'HIST1H2AC', 'NAP1L1']

## plot these 3 first

runlist = main_genes
runlist_name = 'main_genes'

sc.pl.umap(adata_wang, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_no_impute.pdf')) ##
sc.pl.umap(adata_wang_imputed, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_impute.pdf')) ##

runlist = dataset_specific1
runlist_name = 'dataset_specific1'

sc.pl.umap(adata_wang, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_no_impute.pdf')) ##
sc.pl.umap(adata_wang_imputed, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_impute.pdf')) ##

runlist = dataset_specific2
runlist_name = 'dataset_specific2'

sc.pl.umap(adata_wang, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_no_impute.pdf')) ##
sc.pl.umap(adata_wang_imputed, show=False, color=runlist, save=('_wang_full_' + runlist_name + '_impute.pdf')) ##

## then plot the rest -- break into plots with max 12 genes

######

## full list of additional genes -- some may be missing: 
#additional_genes = ['GNLY', 'CCL5', 'NKG7', 'B2M', 'TYROBP', 'IL32', 'SSR4', 'CD74', 'HLA-DRA', 'IL7R', 'RPS12', 'LYZ', 'FTL', 'GZMH', 'CST7', 'MS4A1', 'CD79A', 'CD79B', 'HLA-DRB1', 'CD37', 'BANK1', 'CST3', 'LTB', 'MZB1', 'IGJ', 'FTH1', 'CD14', 'CSF1R', 'IL2RA', 'BATF', 'CX3CR1', 'CXCR3', 'CXCR4', 'CXCR5', 'IRF4', 'HAVCR2', 'NFATC1', 'CCR7', 'SELL', 'TNFRSF9', 'TBX21', 'EOMES', 'PRDM1', 'SLAMF6', 'CD101', 'CCR5', 'IFNG', 'TOX2', 'CD19', 'PTPRC', 'KLRK1', 'TNFRSF17', 'CD80', 'CD86', 'KIT', 'CD34', 'FLT3', 'SIRPA', 'NCAM1', 'FUT4', 'GATA1', 'GATA2', 'GATA3', 'TNF', 'LY6A', 'KLRB1', 'MERTK', 'CCR3', 'TLR7', 'TLR9', 'XCR1', 'CLEC9A', 'CD27', 'CD38', 'ADGRE1', 'FCGR3', 'ZAP70', 'LCK', 'PLCG1', 'TET2', 'LYZ2', 'PEAR1', 'SPI1', 'CEBPA', 'VSIR', 'ITGAM', 'ITGAX', 'ZBTB46', 'CCR9', 'VCAN']

## check for missing genes: 
## [x for x in additional_genes if x not in adata_wang.var.index]
## ['LY6A', 'ADGRE1', 'FCGR3', 'LYZ2', 'VSIR']

## remove these missing genes
additional_genes = ['GNLY', 'CCL5', 'NKG7', 'B2M', 'TYROBP', 'IL32', 'SSR4', 'CD74', 'HLA-DRA', 'IL7R', 'RPS12', 'LYZ', 'FTL', 'GZMH', 'CST7', 'MS4A1', 'CD79A', 'CD79B', 'HLA-DRB1', 'CD37', 'BANK1', 'CST3', 'LTB', 'MZB1', 'IGJ', 'FTH1', 'CD14', 'CSF1R', 'IL2RA', 'BATF', 'CX3CR1', 'CXCR3', 'CXCR4', 'CXCR5', 'IRF4', 'HAVCR2', 'NFATC1', 'CCR7', 'SELL', 'TNFRSF9', 'TBX21', 'EOMES', 'PRDM1', 'SLAMF6', 'CD101', 'CCR5', 'IFNG', 'TOX2', 'CD19', 'PTPRC', 'KLRK1', 'TNFRSF17', 'CD80', 'CD86', 'KIT', 'CD34', 'FLT3', 'SIRPA', 'NCAM1', 'FUT4', 'GATA1', 'GATA2', 'GATA3', 'TNF', 'KLRB1', 'MERTK', 'CCR3', 'TLR7', 'TLR9', 'XCR1', 'CLEC9A', 'CD27', 'CD38', 'ZAP70', 'LCK', 'PLCG1', 'TET2', 'PEAR1', 'SPI1', 'CEBPA', 'ITGAM', 'ITGAX', 'ZBTB46', 'CCR9', 'VCAN']

## set(list(additional_genes)) ## 85 ## 85 / 12 ## 7 full sets then remainder ## set counter to 7

full_sets = len(set(list(additional_genes))) // 12

#####

## break up gene list into parts (8 plots): 

gene_lists = []

counter = 0
endpoint = 0

while counter < full_sets: 
    a = counter * 12
    subset = additional_genes[a:(a+12)]
    endpoint = a+12
    gene_lists.append(subset)
    counter += 1

## if remainder: 
gene_lists.append(additional_genes[endpoint:])

## generate plots

counter2 = 0

for list_x in gene_lists: 

    sc.pl.umap(adata_wang, show=False, color=list_x, save=('_wang_full_' + 'additional_list'  + str(counter2) + '_no_impute.pdf'))
    sc.pl.umap(adata_wang_imputed, show=False, color=list_x, save=('_wang_full_' + 'additional_list' + str(counter2) + '_impute.pdf'))

    counter2 += 1
