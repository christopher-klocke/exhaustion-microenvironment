"""



"""

## import statements 

import scanpy as sc
import numpy as np
import os

## setup 

random1 = 1

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_li_runA'

os.chdir(plots_store)

## file paths 

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

li_just_T_filename = 'li_T_refined.h5ad'

li_T_dpt_filename = 'li_T_dpt2.h5ad'

## load the data

adata_li_T = sc.read_h5ad(li_path + li_just_T_filename)

## process the data

adata_li_T.uns['iroot'] = np.flatnonzero(adata_li_T.obs['leiden'] == '2')[0]

sc.tl.diffmap(adata_li_T, random_state=random1)

sc.tl.dpt(adata_li_T)

## GENERATE PAGA MAP
sc.pl.paga(adata_li_T, show=False, save='_li_t_v1.pdf')

#gene_list = ['TCF7', 'TOX', 'LAG3', 'TIGIT', 'PDCD1']

#sc.pl.diffmap(adata_li_T, color='patient', save='_li_t_patient.pdf', show=False)

#sc.pl.diffmap(adata_li_T, color=gene_list)


## NEED TO SAVE ANNDATA OBJECT WITH DPT ADDED TO LOAD INTO 'set_threshold.py" 

adata_li_T.write_h5ad(li_path + li_T_dpt_filename)
