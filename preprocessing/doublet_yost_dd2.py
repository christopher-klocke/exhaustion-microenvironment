"""
run DoubletDetection method on Yost dataset

https://github.com/JonathanShor/DoubletDetection#running-doubletdetection

run with 'doubletdetection' conda env (instructions in 'roadmap')

raw_counts is a scRNA-seq count matrix (cells by genes), and is array-like
labels is a 1-dimensional numpy ndarray with the value 1 representing a detected doublet, 0 a singlet, and np.nan an ambiguous cell
scores is a 1-dimensional numpy ndarray representing a score for how likely a cell is to be a doublet. The score is used to create the labels

need to correct by sample

save-as from 'doublet_li_dd2.py'
"""

## import statements 

import scanpy as sc
import doubletdetection as dd

## setup 

dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
yost_dir = 'data/yost_data/'  # updated for local
yost_adata_filename = 'yost_viz_ready.h5ad'
yost_doublets_adata_filename = 'yost_doublets_adata.h5ad'
sample_id = 'sample'

## load the data

adata = sc.read_h5ad(dir + yost_dir + yost_adata_filename)

## main function definition
def main(): 
    adata.obs['cell_ids_hold'] = adata.obs.index
    anndata_objects = []
    for p in list(set(adata.obs[sample_id])):
        adata_x = adata[adata.obs[sample_id] == p]
        raw_counts = adata_x.X # raw_counts is a cells by genes count matrix
        clf = dd.BoostClassifier(clustering_algorithm='leiden')
        doublets = clf.fit(raw_counts).predict()
        doublet_score = clf.doublet_score() # higher means more likely to be doublet
        adata_x.obs["doublet"] = doublets
        adata_x.obs["doublet_score"] = doublet_score
        anndata_objects.append(adata_x)

    first_one = anndata_objects.pop(0)
    adata_concat = first_one.concatenate(anndata_objects, join='outer', fill_value=0)
    adata_concat.obs.index = adata_concat.obs['cell_ids_hold']
    del adata_concat.obs['cell_ids_hold']

    adata_concat.write_h5ad(dir + yost_dir + yost_doublets_adata_filename)

## run main function 
if __name__ == "__main__": 
    main()
