"""

viz dd results -- 

from 'doublet_wang_dd2.py' 


"""

## import statements 

import scanpy as sc
import os 

## setup 

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_wang_runA'

os.chdir(plots_store)

wang_dir = '/Users/klockec/Documents/data/wang_data/'

wang_doublets_adata_filename = 'wang_doublets_adata.h5ad'

## load the data

adata = sc.read_h5ad(wang_dir + wang_doublets_adata_filename)

## process the data

sc.pl.umap(adata, color=["doublet", "doublet_score"], show=False, cmap='viridis_r', save='_wang_doublets.pdf')

#sc.pl.umap(adata, color='leiden')

#sc.pl.violin(adata, "doublet_score")

