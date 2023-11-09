"""



"""

## import statements 

import scanpy as sc
import numpy as np
import os 

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/wang_results_plots'

os.chdir(plots_store)

## file paths 

wang_dir = '/Users/klockec/Documents/data/wang_data/'

wang_just_T_filename = 'wang_T_viz_ready.h5ad'

wang_T_dpt_filename = 'wang_T_dpt.h5ad'

## load the data

adata_wang_T = sc.read_h5ad(wang_dir + wang_just_T_filename)

## process the data

adata_wang_T.uns['iroot'] = np.flatnonzero(adata_wang_T.obs['leiden'] == '1')[0]

sc.tl.diffmap(adata_wang_T, random_state=random1)

sc.tl.dpt(adata_wang_T)



gene_list = ['TCF7', 'TOX', 'LAG3', 'TIGIT', 'PDCD1']

sc.pl.diffmap(adata_wang_T, color='batch', save='_wang_t_patient.pdf', show=False)

#sc.pl.diffmap(adata_wang_T, color='dpt_pseudotime')

#sc.pl.diffmap(adata_wang_T, color=gene_list)


## NEED TO SAVE ANNDATA OBJECT WITH DPT ADDED TO LOAD INTO 'set_threshold.py" 

adata_wang_T.write_h5ad(wang_dir + wang_T_dpt_filename)
