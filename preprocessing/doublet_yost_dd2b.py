"""
save-as of 'doublet_yost_dd2.py'

run DoubletDetection method on Yost dataset

https://github.com/JonathanShor/DoubletDetection#running-doubletdetection

run with 'doubletdetection' conda env (instructions in 'roadmap')

raw_counts is a scRNA-seq count matrix (cells by genes), and is array-like
labels is a 1-dimensional numpy ndarray with the value 1 representing a detected doublet, 0 a singlet, and np.nan an ambiguous cell
scores is a 1-dimensional numpy ndarray representing a score for how likely a cell is to be a doublet. The score is used to create the labels

need to correct by sample

input files from 'analysis_process1b_yost.py'
"""

## import statements
import scanpy as sc
import doubletdetection as dd

## setup
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # file path for local
yost_bcc_full_filename = 'data/yost_data/yost_bcc_full_processed.h5ad'
yost_scc_full_filename = 'data/yost_data/yost_scc_full_processed.h5ad'
yost_bcc_doublets_adata_filename = 'data/yost_data/yost_bcc_doublets_adata.h5ad'
yost_scc_doublets_adata_filename = 'data/yost_data/yost_scc_doublets_adata.h5ad'
sample_id = 'sample'

## load the data
adata_bcc = sc.read_h5ad(dir + yost_bcc_full_filename)
adata_scc = sc.read_h5ad(dir + yost_scc_full_filename)

## main function definition
def main():
    ## bcc data
    adata_bcc.obs['cell_ids_hold'] = adata_bcc.obs.index
    anndata_objects_bcc = []
    for p in list(set(adata_bcc.obs[sample_id])):
        adata_x = adata_bcc[adata_bcc.obs[sample_id] == p]
        raw_counts = adata_x.X # raw_counts is a cells by genes count matrix
        clf = dd.BoostClassifier(clustering_algorithm='leiden')
        doublets = clf.fit(raw_counts).predict()
        doublet_score = clf.doublet_score() # higher means more likely to be doublet
        adata_x.obs["doublet"] = doublets
        adata_x.obs["doublet_score"] = doublet_score
        anndata_objects_bcc.append(adata_x)

    first_one_bcc = anndata_objects_bcc.pop(0)
    adata_concat_bcc = first_one_bcc.concatenate(anndata_objects_bcc, join='outer', fill_value=0)
    adata_concat_bcc.obs.index = adata_concat_bcc.obs['cell_ids_hold']
    del adata_concat_bcc.obs['cell_ids_hold']

    adata_concat_bcc.write_h5ad(dir + yost_bcc_doublets_adata_filename)

    ## scc data
    adata_scc.obs['cell_ids_hold'] = adata_scc.obs.index
    anndata_objects_scc = []
    for p in list(set(adata_scc.obs[sample_id])):
        adata_x = adata_scc[adata_scc.obs[sample_id] == p]
        raw_counts = adata_x.X  # raw_counts is a cells by genes count matrix
        clf = dd.BoostClassifier(clustering_algorithm='leiden')
        doublets = clf.fit(raw_counts).predict()
        doublet_score = clf.doublet_score()  # higher means more likely to be doublet
        adata_x.obs["doublet"] = doublets
        adata_x.obs["doublet_score"] = doublet_score
        anndata_objects_scc.append(adata_x)

    first_one_scc = anndata_objects_scc.pop(0)
    adata_concat_scc = first_one_scc.concatenate(anndata_objects_scc, join='outer', fill_value=0)
    adata_concat_scc.obs.index = adata_concat_scc.obs['cell_ids_hold']
    del adata_concat_scc.obs['cell_ids_hold']

    adata_concat_scc.write_h5ad(dir + yost_scc_doublets_adata_filename)

## run main function 
if __name__ == "__main__": 
    main()

