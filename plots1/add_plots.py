"""

generate and save new plots on the fly, as needed

"""

## import statements

import scanpy as sc
import os

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_new'

os.chdir(plots_store)

## file paths

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

li_main_filename = 'li_viz_ready.h5ad'

li_imputed_filename = 'li_imputed.h5ad'

yost_path = '/Users/klockec/Documents/data/yost_data/'

yost_main_filename = 'yost_viz_ready.h5ad'

yost_imputed_filename = 'yost_imputed.h5ad'

wang_path = '/Users/klockec/Documents/data/wang_data/'

wang_main_filename = 'wang_viz_ready.h5ad'

wang_imputed_filename = 'wang_imputed.h5ad'

## load the data

adata_li = sc.read_h5ad(li_path + li_main_filename) ############

#adata_li_imputed = sc.read_h5ad(li_path + li_imputed_filename)

adata_yost = sc.read_h5ad(yost_path + yost_main_filename) ###########

#adata_yost_imputed = sc.read_h5ad(yost_path + yost_imputed_filename)

adata_wang = sc.read_h5ad(wang_path + wang_main_filename) ##########

#adata_wang_imputed = sc.read_h5ad(wang_path + wang_imputed_filename)

######################

## SELECT CORRECT DATASET, IMPUTED OR NOT, GENE(S) TO PLOT, AND OUTPUT FILENAME; THEN RUN SCRIPT AND OPEN PLOT 

gene_list = ['TCF7']

sc.pl.umap(adata_li, show=True, color=gene_list, size=15)
#sc.pl.umap(ANNDATA, show=True, color=gene_list, save=('_DATASET_full_'  + 'FILL' + '_no_impute.pdf'))
