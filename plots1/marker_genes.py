"""

calculate marker genes for each cluster -- for each dataset

"""

## import statements

import scanpy as sc
import os

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img'

os.chdir(plots_store)

## file paths

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

yost_path = '/Users/klockec/Documents/data/yost_data/'

wang_path = '/Users/klockec/Documents/data/wang_data/'

li_filename = 'li_viz_ready.h5ad'

yost_filename = 'yost_viz_ready.h5ad'

wang_filename = 'wang_viz_ready.h5ad'

## load the data

adata_li = sc.read_h5ad(li_path + li_filename)

adata_yost = sc.read_h5ad(yost_path + yost_filename)

adata_wang = sc.read_h5ad(wang_path + wang_filename)

## process the data

adata_li.uns['log1p'] = {'base': None}

adata_yost.uns['log1p'] = {'base': None}

adata_wang.uns['log1p'] = {'base': None}

sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

sc.tl.rank_genes_groups(adata_yost, 'leiden', method='t-test')

sc.tl.rank_genes_groups(adata_wang, 'leiden', method='t-test')

## save plots

sc.pl.rank_genes_groups(adata_li, show=False, save='_li_cluster_marks.pdf')

sc.pl.rank_genes_groups(adata_yost, show=False, save='_yost_cluster_marks.pdf')

sc.pl.rank_genes_groups(adata_wang, show=False, save='_wang_cluster_marks.pdf')
