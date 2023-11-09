"""


"""


## import statements

import scanpy as sc
from cell_separator import separate_score, separate_plot_histograms

## file paths

wang_dir = '/Users/klockec/Documents/data/wang_data/'

wang_just_T_filename = 'wang_T_viz_ready.h5ad'

wang_just_T_imputed_filename = 'wang_T_imputed.h5ad'

## load the data

adata_wang_T = sc.read_h5ad(wang_dir + wang_just_T_imputed_filename)

## process the data 

#genes_dict_cd8_nk = {'higher' : ['CD3E', 'CD3G', 'CD8A', 'CD8B', ], 'lower' : ['GNLY', 'NKG7', 'KLRK1', 'CST7']}  ########
genes_dict_cd8_nk = {'higher' : ['CD3E', 'CD3G', 'CD8A', 'CD8B', ], 'lower' : ['GNLY', 'NKG7', 'KLRK1', 'CST7']}

clusters_list = ['0', '3']

run_scores = separate_score(adata=adata_wang_T, clusters=clusters_list, genes_dict=genes_dict_cd8_nk, save=True, savefile='/Users/klockec/Documents/data/wang_data/split1.txt')

#print(run_scores)

separate_plot_histograms(file='/Users/klockec/Documents/data/wang_data/split1.txt')
