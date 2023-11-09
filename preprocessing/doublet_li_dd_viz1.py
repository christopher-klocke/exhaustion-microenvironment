"""

viz dd results -- 

from 'doublet_li_dd2.py' 


"""

## import statements
import scanpy as sc
import os 

## setup
plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_li_runA'
os.chdir(plots_store)
li_dir = '/Users/klockec/Documents/data/li_data/'
li_doublets_adata_filename = 'li_doublets_adata.h5ad'

## load the data
adata = sc.read_h5ad(li_dir + li_doublets_adata_filename)

## process the data
sc.pl.umap(adata, color=["doublet", "doublet_score"], show=False, cmap='viridis_r', save='_li_doublets.pdf')

#sc.pl.umap(adata, color='leiden')

#sc.pl.violin(adata, "doublet_score")

