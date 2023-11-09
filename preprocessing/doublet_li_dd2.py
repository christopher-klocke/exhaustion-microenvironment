"""
run DoubletDetection method on Li dataset 

https://github.com/JonathanShor/DoubletDetection#running-doubletdetection

run with 'doubletdetection' conda env (instructions in 'roadmap')

raw_counts is a scRNA-seq count matrix (cells by genes), and is array-like
labels is a 1-dimensional numpy ndarray with the value 1 representing a detected doublet, 0 a singlet, and np.nan an ambiguous cell
scores is a 1-dimensional numpy ndarray representing a score for how likely a cell is to be a doublet. The score is used to create the labels

need to correct by sample
"""

## import statements
import scanpy as sc
import doubletdetection as dd

## setup
li_dir = '/Users/klockec/Documents/data/li_data/'
li_adata_filename = 'li_viz_ready.h5ad'
li_doublets_adata_filename = 'li_doublets_adata.h5ad'
#li_dd1_filename = 'li_dd1_save.csv'

## load the data
adata = sc.read_h5ad(li_dir + li_adata_filename)

## main function definition
def main():
    #doublet_dict = {}

    ####### ADD SECTION TO FIX CELL ID ISSUE
    adata.obs['cell_ids_hold'] = adata.obs.index

    anndata_objects = []
    for p in list(set(adata.obs['patient'])): 
        adata_x = adata[adata.obs['patient'] == p]
        raw_counts = adata_x.X # raw_counts is a cells by genes count matrix
        #cell_ids = adata_x.obs.index
        clf = dd.BoostClassifier(clustering_algorithm='leiden')
        doublets = clf.fit(raw_counts).predict()
        doublet_score = clf.doublet_score() # higher means more likely to be doublet
        adata_x.obs["doublet"] = doublets
        adata_x.obs["doublet_score"] = doublet_score
        anndata_objects.append(adata_x)

    first_one = anndata_objects.pop(0)
    adata_concat = first_one.concatenate(anndata_objects, join='outer', fill_value=0)

    ######### ADD SECTION TO FIX CELL ID ISSUE ########
    adata_concat.obs.index = adata_concat.obs['cell_ids_hold']
    del adata_concat.obs['cell_ids_hold']

    adata_concat.write_h5ad(li_dir + li_doublets_adata_filename)

## run main function 
if __name__ == "__main__": 
    main()
